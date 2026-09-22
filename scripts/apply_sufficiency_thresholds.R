#!/usr/bin/env Rscript
# Applies the five data-sufficiency thresholds to output/sufficiency_v1/sufficiency_raw.csv.
#
# Usage (from the repository root):
#   Rscript scripts/apply_sufficiency_thresholds.R
#   Rscript scripts/apply_sufficiency_thresholds.R --coverage-rule b3gbi
#
# Criteria (a combination passes only if it passes all five):
#   density            density_km2         >= 0.25 occurrences / km2
#   Chao2 completeness chao2_completeness  >= 0.70
#   coverage deficit   1 - mean annual coverage <= 0.10
#   temporal           temporal_cha        >= elbow threshold across all combinations
#   precision          precision_ratio     <= 0.10
# A criterion whose value is NA (undefined) counts as a fail.

suppressPackageStartupMessages(library(dplyr))
source("scripts/lib/sufficiency_criteria.R")

arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
IN_FILE  <- arg("--in", "output/sufficiency_v1/sufficiency_raw.csv")
OUT_DIR  <- dirname(IN_FILE)
COVERAGE_RULE <- match.arg(arg("--coverage-rule", "min_cells"), c("min_cells", "b3gbi"))

THRESH <- list(density = 0.25, chao2 = 0.70, coverage_deficit = 0.10, precision = 0.10)

apply_thresholds <- function(raw, thresh = THRESH, coverage_rule = COVERAGE_RULE,
                             temporal_threshold = NULL) {
  x <- raw %>% filter(status == "computed")
  cov_def <- if (coverage_rule == "min_cells") x$coverage_deficit_min_cells else x$coverage_deficit_b3gbi
  if (is.null(temporal_threshold)) temporal_threshold <- elbow_threshold(x$temporal_cha)
  pass <- function(v, ok) !is.na(v) & ok
  x %>% mutate(
    coverage_deficit = cov_def,
    pass_density   = pass(density_km2, density_km2 >= thresh$density),
    pass_chao2     = pass(chao2_completeness, chao2_completeness >= thresh$chao2),
    pass_coverage  = pass(coverage_deficit, coverage_deficit <= thresh$coverage_deficit),
    pass_temporal  = pass(temporal_cha, temporal_cha >= temporal_threshold),
    pass_precision = pass(precision_ratio, precision_ratio <= thresh$precision),
    pass_all = pass_density & pass_chao2 & pass_coverage & pass_temporal & pass_precision,
    temporal_threshold = temporal_threshold,
    coverage_rule = coverage_rule)
}

if (sys.nframe() == 0) {
  raw <- read.csv(IN_FILE, stringsAsFactors = FALSE)
  res <- apply_thresholds(raw)
  cat(sprintf("Coverage rule: %s   temporal (elbow) threshold: %.4f\n",
              COVERAGE_RULE, res$temporal_threshold[1]))
  cat(sprintf("Evaluated combinations: %d (of %d site x taxon rows; the rest had < %d rows)\n",
              nrow(res), nrow(raw), MIN_ROWS))
  for (p in c("pass_density", "pass_chao2", "pass_coverage", "pass_temporal", "pass_precision", "pass_all"))
    cat(sprintf("  %-15s %6d\n", p, sum(res[[p]])))
  cat(sprintf("Passing all five: %d combinations at %d unique sites\n",
              sum(res$pass_all), n_distinct(res$site_id[res$pass_all])))
  print(as.data.frame(res %>% filter(pass_all) %>% count(subset, name = "passing")))
  write.csv(res, file.path(OUT_DIR, "sufficiency_results.csv"), row.names = FALSE)
  write.csv(res %>% filter(pass_all) %>% select(continent, country, site_id, subset),
            file.path(OUT_DIR, "sites_passing.csv"), row.names = FALSE)
}
