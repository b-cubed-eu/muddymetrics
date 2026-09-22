#!/usr/bin/env Rscript
# Trend summaries of the v1 indicators, one row per site x taxon x indicator.
#
# Usage (from the repository root, after run_indicators.R):
#   Rscript scripts/summarise_indicators.R --grid 1km --level 0.9
#   Rscript scripts/summarise_indicators.R --grid 1km --level 0.9 --ci hill_ci_1km_b100.csv
#
# For each combination and indicator: number of years, mean, linear trend
# (slope and p value, as in D6.1: ordinary least squares on annual values,
# direction "up"/"down" when p <= 0.1), and the same trend restricted to the
# MAGV stable window when there is one. Hill estimates that iNEXT could not
# compute at the requested coverage are left out. With --ci, the Hill trends
# are also fitted by weighted least squares, each year weighted by
# 1 / SE^2 with SE taken from the bootstrap interval.
#
# Output: output/indicators_v1/trends_<grid>_<level>.csv and
#         trend_direction_<grid>_<level>.csv

suppressPackageStartupMessages(library(data.table))
arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
GRID <- match.arg(arg("--grid", "1km"), c("100m", "1km"))
LEVEL <- as.numeric(arg("--level", "0.9"))
CI_FILE <- arg("--ci")
OUT <- "output/indicators_v1"

simple <- fread(file.path(OUT, "simple.csv"))[grid == GRID]
hill <- fread(file.path(OUT, "hill.csv"))[grid == GRID & abs(level - LEVEL) < 1e-9 & abs(SC - LEVEL) < 1e-3]
magv <- fread(file.path(OUT, "magv_summary.csv"))

long <- suppressWarnings(rbind(
  melt(simple[, .(site_id, subset, year, total_occ, obs_richness, cum_richness, coverage)],
       id.vars = c("site_id", "subset", "year"), variable.name = "indicator", variable.factor = FALSE),
  hill[, .(site_id, subset, year, indicator = paste0("hill", q), value = qD)]))
if (!is.null(CI_FILE)) {
  ci <- fread(file.path(OUT, CI_FILE))[grid == GRID & abs(level - LEVEL) < 1e-9 & ci_status == "ok" & abs(SC - LEVEL) < 1e-3]
  ci[, se := (qD.UCL - qD.LCL) / (2 * qnorm(0.975))]
  se_tab <- ci[se > 0, .(site_id, subset, year, indicator = paste0("hill", q), se)]
  long <- merge(long, se_tab, by = c("site_id", "subset", "year", "indicator"), all.x = TRUE)
} else long[, se := NA_real_]

fit <- function(v, y, w = NULL) {
  if (length(v) < 3 || length(unique(v)) < 2) return(list(slope = NA_real_, p = NA_real_))
  co <- suppressWarnings(summary(if (is.null(w)) lm(v ~ y) else lm(v ~ y, weights = w)))$coefficients
  if (nrow(co) < 2) return(list(slope = NA_real_, p = NA_real_))
  list(slope = co[2, 1], p = co[2, 4])
}
direction <- function(slope, p) fifelse(is.na(p), NA_character_, fifelse(p <= 0.1, fifelse(slope > 0, "up", "down"), "none"))

long <- merge(long, magv[, .(site_id, subset, stable_run_start, stable_run_end, has_usable_window)], by = c("site_id", "subset"), all.x = TRUE)
tr <- long[order(year), {
  f <- fit(value, year)
  st <- has_usable_window[1] %in% TRUE & !is.na(stable_run_start[1])
  fs <- if (st) { k <- year >= stable_run_start[1] & year <= stable_run_end[1]; fit(value[k], year[k]) } else list(slope = NA_real_, p = NA_real_)
  okw <- !is.na(se)
  fw <- if (sum(okw) >= 3) fit(value[okw], year[okw], 1 / se[okw]^2) else list(slope = NA_real_, p = NA_real_)
  .(n_years = .N, first_year = min(year), last_year = max(year), mean_value = mean(value),
    slope = f$slope, p = f$p, slope_stable = fs$slope, p_stable = fs$p,
    n_years_with_se = sum(okw), slope_weighted = fw$slope, p_weighted = fw$p)
}, by = .(site_id, subset, indicator)]
tr[, `:=`(direction = direction(slope, p), direction_weighted = direction(slope_weighted, p_weighted))]
tag <- paste0(GRID, "_", LEVEL)
fwrite(tr, file.path(OUT, paste0("trends_", tag, ".csv")))
dir <- tr[!is.na(direction), .(combinations = .N, up = sum(direction == "up"), none = sum(direction == "none"),
                               down = sum(direction == "down")), by = .(subset, indicator)]
fwrite(dir, file.path(OUT, paste0("trend_direction_", tag, ".csv")))
cat(sprintf("Grid %s, coverage %.2f: %d trend rows written to %s\n", GRID, LEVEL, nrow(tr), file.path(OUT, paste0("trends_", tag, ".csv"))))
