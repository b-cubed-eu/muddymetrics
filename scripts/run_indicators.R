#!/usr/bin/env Rscript
# Runs worker_indicators.R over every site, for every site x taxon combination
# the sufficiency screen evaluated, and merges the results.
#
# Usage (from the repository root, after run_sufficiency.R and
# apply_sufficiency_thresholds.R):
#   Rscript scripts/run_indicators.R --cores 8                         # point estimates, all combinations
#   Rscript scripts/run_indicators.R --cores 8 --mode ci --nboot 100 --grids 1km   # intervals, passing combinations
#   Rscript scripts/run_indicators.R --screen output/sufficiency_v1b_1km/sufficiency_raw.csv \
#     --coverage-rule informative --out output/indicators_v1b --sites passing_sites.txt
#   Rscript scripts/run_indicators.R --mode ci --nboot 1000 --grids 1km --years 2024-2025 --out output/indicators_y2425
#
# --screen and --coverage-rule set which screen results decide the combinations
# (default: the 100 m screen with the min_cells rule); --out the output folder.
# --combos: a CSV with site_id and subset columns; in mode "ci", only these
# combinations get intervals (for adding combinations to an earlier run).
# --years: mode "ci" only, e.g. 2024-2025: intervals for these years only (for
# adding years to an earlier run; each year is estimated on its own).
# --split-hours H: mode "ci" only; splits each combination by year into parts of
# about H hours, run as separate processes, so that one large combination does not
# hold up the run (each year and order has its own seed, so the merged result is
# the same as one run over all years). Part sizes come from the seconds of an
# earlier run given with --prior <hill_ci file>; combinations without an earlier
# time get a time predicted from their size in the screen (or, if that is not
# possible, --split-default parts [8]). Parts go to <out>/parts/ and are
# merged into the usual per-site files when all parts of a site are done.
#   Rscript scripts/run_indicators.R --mode ci --nboot 1000 --grids 1km --split-hours 0.5 \
#     --prior output/indicators_v1_b1000/hill_ci_1km_b1000.csv ...
#
# Mode "point": every combination the screen evaluated gets point estimates.
# Mode "ci": combinations that pass all five criteria under any of the three
# precision rules (mean, median, D6.1) get bootstrap intervals for the Hill
# numbers.
#
# Output in output/indicators_v1/ (or --out): simple.csv, hill.csv, hill_ci_<grids>_b<nboot>.csv,
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
YEARS <- arg("--years", "")
if (nzchar(YEARS) && mode != "ci") stop("--years can only be used with --mode ci")
ci_suffix <- paste0("_hill_ci_", gsub(",", "_", grids), "_b", nboot, if (nzchar(YEARS)) paste0("_y", YEARS) else "", ".csv")
sites_file <- arg("--sites")
OUT <- arg("--out", "output/indicators_v1")
SCREEN <- arg("--screen", "output/sufficiency_v1/sufficiency_raw.csv")
COV_RULE <- arg("--coverage-rule", "min_cells")
COMBOS <- arg("--combos")
SPLIT_H <- as.numeric(arg("--split-hours", "0"))
PRIOR <- arg("--prior")
SPLIT_DEFAULT <- as.integer(arg("--split-default", "8"))
if (SPLIT_H > 0 && mode != "ci") stop("--split-hours can only be used with --mode ci")
if (SPLIT_H > 0 && nzchar(YEARS)) stop("--split-hours cannot be combined with --years")

raw <- read.csv(SCREEN, stringsAsFactors = FALSE)
evaluated <- raw[raw$status == "computed", c("site_id", "subset")]
pass <- rbindlist(lapply(c("mean", "median", "d61"), function(pr) {
  r <- apply_thresholds(raw, coverage_rule = COV_RULE, precision_rule = pr)
  data.table(site_id = r$site_id[r$pass_all], subset = r$subset[r$pass_all], precision_rule = pr)
}))
ci_combos <- unique(pass[, .(site_id, subset)])
if (!is.null(COMBOS)) ci_combos <- unique(fread(COMBOS)[, .(site_id, subset)])

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

if (SPLIT_H > 0) {
  source("scripts/lib/indicators_v1.R")   # ci_job_plan(), merge_ci_parts()
  dir.create(file.path(OUT, "parts"), recursive = TRUE, showWarnings = FALSE)
  combos_todo <- ci_combos[site_id %in% todo$site_id]
  plan <- ci_job_plan(combos_todo, prior = if (!is.null(PRIOR)) fread(PRIOR) else NULL, hours = SPLIT_H, default_k = SPLIT_DEFAULT,
                      size = raw[raw$status == "computed", ])
  plan[, path := todo$path[match(site_id, todo$site_id)]]
  plan[, prefix := file.path(OUT, "parts", paste0(site_id, "_", subset))]
  plan[, file := paste0(prefix, sub("\\.csv$", "", ci_suffix), "_y", years, ".csv")]
  fwrite(plan[, .(site_id, subset, k, part, years, est_hours, predicted, est_part_hours)], file.path(OUT, "split_plan.csv"))
  jt <- plan[!file.exists(file)]
  cat(sprintf("Split into %d parts (%d already done); expected longest part %.1f h, total %.0f h of computing.\n",
              nrow(plan), nrow(plan) - nrow(jt), max(plan$est_part_hours), sum(plan$est_part_hours)))
  if (nrow(jt)) {
    cl <- makePSOCKcluster(n_cores)
    st <- parLapplyLB(cl, seq_len(nrow(jt)), function(i, jt, OUT, nboot, grids) {
      j <- jt[i, ]
      log <- file.path(OUT, "logs", paste0(j$site_id, "_", j$subset, "_", j$years, ".log"))
      code <- system2("Rscript", c("scripts/worker_indicators.R", j$site_id, shQuote(j$path), shQuote(j$prefix),
                                   j$subset, "ci", nboot, grids, j$years), stdout = log, stderr = log)
      ok <- code == 0 && file.exists(j$file)
      if (ok) unlink(log)
      data.frame(site_id = j$site_id, subset = j$subset, years = j$years, ok = ok, exit_code = code)
    }, jt = as.data.frame(jt), OUT = OUT, nboot = nboot, grids = grids, chunk.size = 1)
    stopCluster(cl)
    st <- rbindlist(st)
    cat(sprintf("Finished. %d parts succeeded, %d failed.\n", sum(st$ok), sum(!st$ok)))
    if (any(!st$ok)) fwrite(st[!st$ok], file.path(OUT, "failures_ci_parts.csv"))
  }
  # sites whose parts are all done get the per-site file an unsplit run would have written
  for (s in unique(plan$site_id)) {
    p <- plan[site_id == s]
    if (!all(file.exists(p$file))) next
    dest <- paste0(file.path(OUT, "sites", s), ci_suffix)
    m <- merge_ci_parts(p$file, subset_order = strsplit(todo$subsets[todo$site_id == s], ",")[[1]])
    fwrite(m, paste0(dest, ".tmp")); file.rename(paste0(dest, ".tmp"), dest)
  }
  cat(sprintf("Sites complete: %d of %d.\n", sum(file.exists(paste0(file.path(OUT, "sites", unique(plan$site_id)), ci_suffix))), uniqueN(plan$site_id)))
} else if (nrow(todo) > 0) {
  cl <- makePSOCKcluster(n_cores)
  st <- parLapplyLB(cl, seq_len(nrow(todo)), function(i, todo, OUT, nboot, mode, grids, ci_suffix, YEARS) {
    j <- todo[i, ]
    log <- file.path(OUT, "logs", paste0(j$site_id, "_", mode, ".log"))
    code <- system2("Rscript", c("scripts/worker_indicators.R", j$site_id, shQuote(j$path), shQuote(j$prefix),
                                 j$subsets, mode, nboot, grids, YEARS), stdout = log, stderr = log)
    ok <- code == 0 && file.exists(paste0(j$prefix, if (mode == "point") "_magv.csv" else ci_suffix))
    if (ok) unlink(log)
    data.frame(site_id = j$site_id, ok = ok, exit_code = code)
  }, todo = as.data.frame(todo), OUT = OUT, nboot = nboot, mode = mode, grids = grids, ci_suffix = ci_suffix, YEARS = YEARS, chunk.size = 1)
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
