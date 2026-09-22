#!/usr/bin/env Rscript
# Runs worker_indicators.R over every site, for every site x taxon combination
# the sufficiency screen evaluated, and merges the results.
#
# Usage (from the repository root, after run_sufficiency.R and
# apply_sufficiency_thresholds.R):
#   Rscript scripts/run_indicators.R --cores 8                         # point estimates, all combinations
#   Rscript scripts/run_indicators.R --cores 8 --mode ci --nboot 100 --grids 1km   # intervals, passing combinations
#
# Mode "point": every combination the screen evaluated gets point estimates.
# Mode "ci": combinations that pass all five criteria under any of the three
# precision rules (mean, median, D6.1) get bootstrap intervals for the Hill
# numbers.
#
# Output in output/indicators_v1/: simple.csv, hill.csv, hill_ci_<grids>_b<nboot>.csv,
# magv_by_year.csv, magv_summary.csv, ci_combinations.csv, per-site files in sites/.
# Resumable: sites with finished files for the chosen mode are skipped.

suppressPackageStartupMessages({ library(parallel); library(data.table) })
source("scripts/lib/sufficiency_criteria.R")
source("scripts/lib/magv.R")
source("scripts/apply_sufficiency_thresholds.R")   # apply_thresholds(); main block does not run when sourced

arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
n_cores <- as.integer(arg("--cores", "8"))
nboot <- as.integer(arg("--nboot", "200"))
mode <- match.arg(arg("--mode", "point"), c("point", "ci"))
grids <- arg("--grids", "100m,1km")
ci_suffix <- paste0("_hill_ci_", gsub(",", "_", grids), "_b", nboot, ".csv")
sites_file <- arg("--sites")
OUT <- "output/indicators_v1"

raw <- read.csv("output/sufficiency_v1/sufficiency_raw.csv", stringsAsFactors = FALSE)
evaluated <- raw[raw$status == "computed", c("site_id", "subset")]
pass <- rbindlist(lapply(c("mean", "median", "d61"), function(pr) {
  r <- apply_thresholds(raw, coverage_rule = "min_cells", precision_rule = pr)
  data.table(site_id = r$site_id[r$pass_all], subset = r$subset[r$pass_all], precision_rule = pr)
}))
ci_combos <- unique(pass[, .(site_id, subset)])

cubes <- discover_site_cubes()
sites <- data.table(site_id = unique(evaluated$site_id))
if (!is.null(sites_file)) sites <- sites[site_id %in% trimws(readLines(sites_file))]
sites$path <- cubes$file_path[match(sites$site_id, cubes$site_id)]
sites$size <- file.info(sites$path)$size
sites$subsets <- vapply(sites$site_id, function(s) paste(evaluated$subset[evaluated$site_id == s], collapse = ","), "")
sites$ci <- vapply(sites$site_id, function(s) paste(ci_combos$subset[ci_combos$site_id == s], collapse = ","), "")
if (mode == "ci") { sites <- sites[ci != ""]; sites$subsets <- sites$ci }
dir.create(file.path(OUT, "sites"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT, "logs"), recursive = TRUE, showWarnings = FALSE)
sites$prefix <- file.path(OUT, "sites", sites$site_id)
done_file <- paste0(sites$prefix, if (mode == "point") "_magv.csv" else ci_suffix)
todo <- sites[!file.exists(done_file)][order(-size)]
cat(sprintf("Mode %s: %d sites (%d evaluated combinations; %d passing under some precision rule%s). %d already done, %d to run on %d cores.\n",
            mode, nrow(sites), nrow(evaluated), nrow(ci_combos), if (mode == "ci") sprintf(", nboot = %d", nboot) else "",
            nrow(sites) - nrow(todo), nrow(todo), n_cores))

if (nrow(todo) > 0) {
  cl <- makePSOCKcluster(n_cores)
  st <- parLapplyLB(cl, seq_len(nrow(todo)), function(i, todo, OUT, nboot, mode, grids, ci_suffix) {
    j <- todo[i, ]
    log <- file.path(OUT, "logs", paste0(j$site_id, "_", mode, ".log"))
    code <- system2("Rscript", c("scripts/worker_indicators.R", j$site_id, shQuote(j$path), shQuote(j$prefix),
                                 j$subsets, mode, nboot, grids), stdout = log, stderr = log)
    ok <- code == 0 && file.exists(paste0(j$prefix, if (mode == "point") "_magv.csv" else ci_suffix))
    if (ok) unlink(log)
    data.frame(site_id = j$site_id, ok = ok, exit_code = code)
  }, todo = as.data.frame(todo), OUT = OUT, nboot = nboot, mode = mode, grids = grids, ci_suffix = ci_suffix, chunk.size = 1)
  stopCluster(cl)
  st <- rbindlist(st)
  cat(sprintf("Finished. %d succeeded, %d failed.\n", sum(st$ok), sum(!st$ok)))
  if (any(!st$ok)) fwrite(st[!st$ok], file.path(OUT, paste0("failures_", mode, ".csv")))
}

merge_part <- function(suffix) {
  f <- paste0(sites$prefix, suffix); f <- f[file.exists(f)]
  rbindlist(lapply(f, fread), fill = TRUE)
}
if (mode == "point") {
  fwrite(merge_part("_simple.csv"), file.path(OUT, "simple.csv"))
  fwrite(merge_part("_hill.csv"), file.path(OUT, "hill.csv"))
  mg <- merge_part("_magv.csv")
  fwrite(mg, file.path(OUT, "magv_by_year.csv"))
  if (nrow(mg)) fwrite(mg[, magv_summary(as.data.frame(.SD)), by = .(site_id, subset)], file.path(OUT, "magv_summary.csv"))
} else fwrite(merge_part(ci_suffix), file.path(OUT, sub("^_", "", ci_suffix)))
fwrite(pass, file.path(OUT, "ci_combinations.csv"))
cat("Merged results in", OUT, "\n")
