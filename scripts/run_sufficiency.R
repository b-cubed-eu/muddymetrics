#!/usr/bin/env Rscript
# Runs worker_sufficiency.R over every site cube in parallel and merges the
# per-site results into output/sufficiency_v1/sufficiency_raw.csv.
#
# Usage (from the repository root):
#   Rscript scripts/run_sufficiency.R                 # all sites, 8 cores
#   Rscript scripts/run_sufficiency.R --cores 4
#   Rscript scripts/run_sufficiency.R --sites pilot_sites.txt   # one site_id per line
#   Rscript scripts/run_sufficiency.R --grid 1km      # screen on 1 km cells
#                                    # (results in output/sufficiency_v1_1km/)
#   Rscript scripts/run_sufficiency.R --out output/my_run   # another output folder
#   Rscript scripts/run_sufficiency.R --last-year 2023     # leave out records after 2023
#
# Resumable: sites that already have a result file are skipped. Each site is
# run in its own R process, so one failing site cannot stop the run; failures
# are listed in failures.csv with the error output kept in logs/.

suppressPackageStartupMessages({ library(dplyr); library(parallel) })
source("scripts/lib/sufficiency_criteria.R")

arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
n_cores <- as.integer(arg("--cores", "8"))
sites_file <- arg("--sites")
grid <- match.arg(arg("--grid", "100m"), c("100m", "1km"))
OUT_DIR <- arg("--out", if (grid == "100m") "output/sufficiency_v1" else "output/sufficiency_v1_1km")
last_year <- arg("--last-year", "Inf")

# Site discovery (scripts/lib/sufficiency_criteria.R), identical to
# run_subgroup_data_sufficiency_parallel.R so the site set is comparable with
# D6.1: one cube per site_id, the largest file where a site appears more than once.
sites <- discover_site_cubes()
cat(sprintf("Discovered %d sites.\n", nrow(sites)))
if (!is.null(sites_file)) {
  wanted <- trimws(readLines(sites_file)); wanted <- wanted[wanted != ""]
  missing <- setdiff(wanted, sites$site_id)
  if (length(missing)) cat("Not found, skipped:", paste(missing, collapse = " "), "\n")
  sites <- sites[sites$site_id %in% wanted, ]
  cat(sprintf("Restricted to %d sites from %s.\n", nrow(sites), sites_file))
}

site_dir <- file.path(OUT_DIR, "sites"); log_dir <- file.path(OUT_DIR, "logs")
dir.create(site_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
sites$out_file <- file.path(site_dir, paste0(sites$site_id, ".csv"))
todo <- sites[!file.exists(sites$out_file), ]
cat(sprintf("%d already done, %d to run on %d cores.\n", nrow(sites) - nrow(todo), nrow(todo), n_cores))

if (nrow(todo) > 0) {
  cl <- makePSOCKcluster(n_cores)
  status <- parLapply(cl, seq_len(nrow(todo)), function(i, todo, log_dir, grid, last_year) {
    s <- todo[i, ]
    log <- file.path(log_dir, paste0(s$site_id, ".log"))
    code <- system2("Rscript", c("scripts/worker_sufficiency.R", shQuote(s$cont), shQuote(s$country),
                                 shQuote(s$site_id), shQuote(s$file_path), shQuote(s$out_file), grid, last_year),
                    stdout = log, stderr = log)
    ok <- code == 0 && file.exists(s$out_file)
    if (ok) unlink(log)   # keep logs only for failures
    data.frame(site_id = s$site_id, ok = ok, exit_code = code)
  }, todo = todo, log_dir = log_dir, grid = grid, last_year = last_year)
  stopCluster(cl)
  status <- bind_rows(status)
  failed <- status[!status$ok, ]
  cat(sprintf("Finished. %d succeeded, %d failed.\n", sum(status$ok), nrow(failed)))
  if (nrow(failed)) write.csv(failed, file.path(OUT_DIR, "failures.csv"), row.names = FALSE)
}

done <- sites$out_file[file.exists(sites$out_file)]
raw <- bind_rows(lapply(done, function(f) read.csv(f, stringsAsFactors = FALSE)))
write.csv(raw, file.path(OUT_DIR, "sufficiency_raw.csv"), row.names = FALSE)
cat(sprintf("Merged %d sites (%d site x taxon rows, %d computed) into %s\n",
            length(done), nrow(raw), sum(raw$status == "computed"),
            file.path(OUT_DIR, "sufficiency_raw.csv")))
