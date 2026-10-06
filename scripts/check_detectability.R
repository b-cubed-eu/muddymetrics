#!/usr/bin/env Rscript
# Smallest detectable trend for each passing site x taxon combination: the
# steady linear change over a window of years that an ordinary least squares
# trend test would detect with the given power at the given significance level,
# expressed as a percentage of the combination's mean value per decade.
#
#   MDT = (t[1 - alpha/2, n-2] + t[power, n-2]) * SE(slope) * 10 / mean * 100
#
# SE(slope) is computed two ways: from the residual scatter around the fitted
# line (observed year-to-year variation), and from the bootstrap standard
# errors of the yearly estimates alone (sampling error only).
#
# Usage (from the repository root, after run_indicators.R in both modes):
#   Rscript scripts/check_detectability.R
#   Rscript scripts/check_detectability.R --grid 1km --level 0.9 --from 2016 --to 2025 \
#     --min-years 8 --alpha 0.05 --power 0.8 --precision-rule median --ci hill_ci_1km_b100.csv
#   Rscript scripts/check_detectability.R --in output/indicators_v1b --combos output/sufficiency_v1b_1km/sites_passing.csv
#
# Hill estimates use the same years as summarise_indicators.R: those reaching the
# coverage level within --max-extrapolation (default 2) times the observed cells,
# in informative years (--hill-years informative, default; or all).
# --in is the indicator folder; --combos a CSV with site_id and subset columns.
#
# Output: output/v1_checks/mdt_per_combo_<rule>.csv and a summary on screen.

suppressPackageStartupMessages(library(data.table))
arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
GRID <- match.arg(arg("--grid", "1km"), c("100m", "1km"))
LEVEL <- as.numeric(arg("--level", "0.9"))
FROM <- as.integer(arg("--from", "2016")); TO <- as.integer(arg("--to", "2025"))
MIN_YEARS <- as.integer(arg("--min-years", "8"))
ALPHA <- as.numeric(arg("--alpha", "0.05")); POWER <- as.numeric(arg("--power", "0.8"))
RULE <- match.arg(arg("--precision-rule", "mean"), c("mean", "median", "d61"))
CI_FILE <- arg("--ci", sprintf("hill_ci_%s_b100.csv", GRID))
BENCHMARK <- 30   # % per decade: a 30% decline in 10 years (IUCN Red List criterion A2, Vulnerable)
IN <- arg("--in", "output/indicators_v1"); OUT <- "output/v1_checks"
COMBOS <- arg("--combos")
HILL_YEARS <- match.arg(arg("--hill-years", "informative"), c("informative", "all"))
MAX_EXTRAP <- as.numeric(arg("--max-extrapolation", "2"))
source("scripts/lib/sufficiency_criteria.R")   # INFORMATIVE_MIN_CELLS, INFORMATIVE_MIN_INCIDENCES
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

combos <- if (is.null(COMBOS)) unique(fread(file.path(IN, "ci_combinations.csv"))[precision_rule == RULE, .(site_id, subset)]) else
  unique(fread(COMBOS)[, .(site_id, subset)])
h <- fread(file.path(IN, "hill.csv"))[grid == GRID & abs(level - LEVEL) < 1e-9 & abs(SC - LEVEL) < 1e-3]
ci <- if (file.exists(file.path(IN, CI_FILE))) {
  x <- fread(file.path(IN, CI_FILE))[grid == GRID & abs(level - LEVEL) < 1e-9 & ci_status == "ok" & abs(SC - LEVEL) < 1e-3]
  x[, se := (qD.UCL - qD.LCL) / (2 * qnorm(0.975))]
} else data.table(site_id = character(), subset = character(), year = integer(), q = integer(), se = numeric())
s <- fread(file.path(IN, "simple.csv"))[grid == GRID]
h <- h[t_over_T <= MAX_EXTRAP]
if (HILL_YEARS == "informative") {
  if (!"n_incidences" %in% names(s)) stop("simple.csv has no n_incidences column: re-run run_indicators.R, or use --hill-years all")
  h <- merge(h, s[n_cells >= INFORMATIVE_MIN_CELLS & n_incidences >= INFORMATIVE_MIN_INCIDENCES, .(site_id, subset, year)],
             by = c("site_id", "subset", "year"))
}

series <- rbind(
  merge(h[q %in% c(0, 1), .(site_id, subset, year, q, v = qD)], ci[, .(site_id, subset, year, q, se)],
        by = c("site_id", "subset", "year", "q"), all.x = TRUE)[, indicator := paste0("hill", q)][, q := NULL],
  s[, .(site_id, subset, year, v = obs_richness, se = NA_real_, indicator = "obs_richness")])
series <- merge(series, combos, by = c("site_id", "subset"))[year >= FROM & year <= TO]

res <- series[order(year), if (.N >= MIN_YEARS && mean(v) > 0) {
  t <- year - mean(year); sxx <- sum(t^2); n <- .N
  k <- qt(1 - ALPHA / 2, n - 2) + qt(POWER, n - 2)
  se_obs <- summary(lm(v ~ year))$sigma / sqrt(sxx)
  se_samp <- if (all(!is.na(se))) sqrt(sum(t^2 * se^2)) / sxx else NA_real_
  .(n_years = n, mean = mean(v), multiplier = k,
    mdt_observed = 100 * k * se_obs * 10 / mean(v), mdt_sampling_only = 100 * k * se_samp * 10 / mean(v))
}, by = .(site_id, subset, indicator)]

tag <- if (is.null(COMBOS)) RULE else basename(dirname(COMBOS))
fwrite(res, file.path(OUT, sprintf("mdt_per_combo_%s.csv", tag)))
cat(sprintf("Smallest detectable linear change, %d-%d, %s grid, %.0f%% coverage, alpha %.2f, power %.2f, >= %d years, combinations: %s\n",
            FROM, TO, GRID, 100 * LEVEL, ALPHA, POWER, MIN_YEARS, if (is.null(COMBOS)) paste("precision rule", RULE) else COMBOS))
print(res[, .(combinations = .N, median_mdt = round(median(mdt_observed)),
              share_can_detect_benchmark = round(mean(mdt_observed <= BENCHMARK), 2),
              median_mdt_sampling_only = round(median(mdt_sampling_only, na.rm = TRUE))), by = indicator])
