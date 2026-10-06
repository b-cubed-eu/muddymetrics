#!/usr/bin/env Rscript
# Applies the five data-sufficiency thresholds to output/sufficiency_v1/sufficiency_raw.csv.
#
# Usage (from the repository root):
#   Rscript scripts/apply_sufficiency_thresholds.R
#   Rscript scripts/apply_sufficiency_thresholds.R --coverage-rule b3gbi
#   Rscript scripts/apply_sufficiency_thresholds.R --coverage-rule informative --min-informative-years 10
#   Rscript scripts/apply_sufficiency_thresholds.R --precision-rule median
#
# Criteria (a combination passes only if it passes all five):
#   density            density_km2         >= 0.25 occurrences / km2
#   Chao2 completeness chao2_completeness  >= 0.70
#   coverage deficit   1 - mean annual coverage <= 0.10
#                      (--coverage-rule min_cells: years with >= 2 occupied cells, default;
#                       b3gbi: b3gbi's matrix-size rule; informative: years with
#                       >= INFORMATIVE_MIN_CELLS cells and >= INFORMATIVE_MIN_INCIDENCES incidences)
#   temporal           temporal_cha        >= elbow threshold across all combinations
#   precision          precision_ratio     <= 0.10
#                      (--precision-rule mean: occurrence-weighted mean uncertainty area, default;
#                       median: occurrence-weighted median; d61: the D6.1 site-level definition)
# With --coverage-rule informative, the coverage criterion also needs at least
# --min-informative-years informative years (default 0: no minimum).
# A criterion whose value is NA (undefined) counts as a fail.

suppressPackageStartupMessages(library(dplyr))
source("scripts/lib/sufficiency_criteria.R")

arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
IN_FILE  <- arg("--in", "output/sufficiency_v1/sufficiency_raw.csv")
OUT_DIR  <- dirname(IN_FILE)
COVERAGE_RULE <- match.arg(arg("--coverage-rule", "min_cells"), c("min_cells", "b3gbi", "informative"))
PRECISION_RULE <- match.arg(arg("--precision-rule", "mean"), c("mean", "median", "d61"))
MIN_INF_YEARS <- as.integer(arg("--min-informative-years", "0"))
PRECISION_COL <- c(mean = "precision_ratio", median = "precision_ratio_median", d61 = "precision_ratio_d61")

THRESH <- list(density = 0.25, chao2 = 0.70, coverage_deficit = 0.10, precision = 0.10)

apply_thresholds <- function(raw, thresh = THRESH, coverage_rule = COVERAGE_RULE,
                             precision_rule = PRECISION_RULE, temporal_threshold = NULL,
                             min_informative_years = MIN_INF_YEARS) {
  x <- raw %>% filter(status == "computed")
  cov_def <- x[[paste0("coverage_deficit_", coverage_rule)]]
  if (is.null(cov_def)) stop("No coverage_deficit_", coverage_rule, " column: re-run the screen")
  if (coverage_rule == "informative") cov_def[x$n_years_informative < min_informative_years] <- NA_real_
  prec <- x[[PRECISION_COL[[precision_rule]]]]
  if (is.null(temporal_threshold)) temporal_threshold <- elbow_threshold(x$temporal_cha)
  pass <- function(v, ok) !is.na(v) & ok
  x %>% mutate(
    coverage_deficit = cov_def,
    precision_used = prec,
    pass_density   = pass(density_km2, density_km2 >= thresh$density),
    pass_chao2     = pass(chao2_completeness, chao2_completeness >= thresh$chao2),
    pass_coverage  = pass(coverage_deficit, coverage_deficit <= thresh$coverage_deficit),
    pass_temporal  = pass(temporal_cha, temporal_cha >= temporal_threshold),
    pass_precision = pass(precision_used, precision_used <= thresh$precision),
    pass_all = pass_density & pass_chao2 & pass_coverage & pass_temporal & pass_precision,
    temporal_threshold = temporal_threshold,
    coverage_rule = coverage_rule,
    precision_rule = precision_rule)
}

if (sys.nframe() == 0) {
  raw <- read.csv(IN_FILE, stringsAsFactors = FALSE)
  res <- apply_thresholds(raw)
  cat(sprintf("Coverage rule: %s   precision rule: %s   temporal (elbow) threshold: %.4f\n",
              COVERAGE_RULE, PRECISION_RULE, res$temporal_threshold[1]))
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
  cat("\nAll rule combinations (combinations passing all five / unique sites):\n")
  for (cr in intersect(c("min_cells", "b3gbi", "informative"), sub("^coverage_deficit_", "", names(raw)))) for (pr in names(PRECISION_COL)) {
    r <- apply_thresholds(raw, coverage_rule = cr, precision_rule = pr)
    cat(sprintf("  coverage %-9s precision %-6s %5d / %d\n", cr, pr, sum(r$pass_all),
                n_distinct(r$site_id[r$pass_all])))
  }
}
