#!/usr/bin/env Rscript
# Runs worker_indicators.R for every site x taxon combination that passed the
# data-sufficiency screen, then merges the MAGV and status tables.
#
# Usage (from the repository root, after apply_sufficiency_thresholds.R):
#   Rscript scripts/run_indicators.R
#   Rscript scripts/run_indicators.R --cores 2 --in output/sufficiency_v1/sites_passing.csv
#
# Output: output/indicators_v1/<site_id>/<subset>_<indicator>_ts.rds, plus
#   output/indicators_v1/magv_by_year.csv, magv_summary.csv and status.csv.
# Resumable: combinations with a status file are skipped. Taxonomic
# distinctness queries the GBIF API (taxize), so the run needs internet access.

suppressPackageStartupMessages({ library(parallel); library(data.table) })
source("scripts/lib/sufficiency_criteria.R")

arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
n_cores <- as.integer(arg("--cores", "2"))
in_file <- arg("--in", "output/sufficiency_v1/sites_passing.csv")
OUT_DIR <- arg("--out", "output/indicators_v1")

if (packageVersion("b3gbi") < "1.0.0") stop("b3gbi >= 1.0.0 is required; found ", packageVersion("b3gbi"))

jobs <- read.csv(in_file, stringsAsFactors = FALSE)
cubes <- discover_site_cubes()
jobs$file_path <- cubes$file_path[match(jobs$site_id, cubes$site_id)]
if (anyNA(jobs$file_path)) stop("No cube found for: ", paste(jobs$site_id[is.na(jobs$file_path)], collapse = " "))
jobs$out_dir <- file.path(OUT_DIR, jobs$site_id)
jobs$status_file <- file.path(jobs$out_dir, paste0(jobs$subset, "_status.csv"))
todo <- jobs[!file.exists(jobs$status_file), ]
cat(sprintf("%d combinations, %d already done, %d to run on %d cores.\n",
            nrow(jobs), nrow(jobs) - nrow(todo), nrow(todo), n_cores))

log_dir <- file.path(OUT_DIR, "logs"); dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
if (nrow(todo) > 0) {
  cl <- makePSOCKcluster(n_cores)
  res <- parLapply(cl, seq_len(nrow(todo)), function(i, todo, log_dir) {
    j <- todo[i, ]
    log <- file.path(log_dir, paste0(j$site_id, "_", j$subset, ".log"))
    code <- system2("Rscript", c("scripts/worker_indicators.R", shQuote(j$continent), shQuote(j$country),
                                 shQuote(j$site_id), shQuote(j$subset), shQuote(j$file_path),
                                 shQuote(j$out_dir)), stdout = log, stderr = log)
    ok <- code == 0 && file.exists(j$status_file)
    if (ok) unlink(log)
    data.frame(site_id = j$site_id, subset = j$subset, ok = ok, exit_code = code)
  }, todo = todo, log_dir = log_dir)
  stopCluster(cl)
  res <- rbindlist(res)
  cat(sprintf("Finished. %d succeeded, %d failed.\n", sum(res$ok), sum(!res$ok)))
  if (any(!res$ok)) fwrite(res[!res$ok], file.path(OUT_DIR, "failures.csv"))
}

merge_files <- function(pattern) {
  f <- list.files(OUT_DIR, pattern = pattern, recursive = TRUE, full.names = TRUE)
  rbindlist(lapply(f, fread), fill = TRUE)
}
fwrite(merge_files("_magv_by_year\\.csv$"), file.path(OUT_DIR, "magv_by_year.csv"))
fwrite(merge_files("_magv_summary\\.csv$"), file.path(OUT_DIR, "magv_summary.csv"))
st <- merge_files("_status\\.csv$")
fwrite(st, file.path(OUT_DIR, "status.csv"))
if (nrow(st)) {
  cat("Indicator runs by CI outcome:\n"); print(table(st$indicator, st$ci))
}
