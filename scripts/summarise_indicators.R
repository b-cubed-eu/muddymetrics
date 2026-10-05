#!/usr/bin/env Rscript
# Trend summaries of the v1 indicators, one row per site x taxon x indicator.
#
# Usage (from the repository root, after run_indicators.R):
#   Rscript scripts/summarise_indicators.R --grid 1km --level 0.9
#   Rscript scripts/summarise_indicators.R --grid 1km --level 0.9 --ci hill_ci_1km_b100.csv
#   Rscript scripts/summarise_indicators.R --grid 1km --level 0.9 --alpha 0.1   # D6.1's level
#
# For each combination and indicator: number of years, mean, and linear trends
# fitted three ways:
#   ols  ordinary least squares on annual values (as in D6.1), also restricted
#        to the MAGV stable window when there is one;
#   gls  generalised least squares with AR(1) errors, fitted by REML
#        (nlme::gls with corAR1 on year);
#   with --ci, the Hill trends are also fitted with each year weighted by
#        1 / SE^2 (SE from the bootstrap interval): wls (least squares) and
#        gls_w (AR(1) errors, variance proportional to SE^2).
# Direction is "up"/"down" when p <= alpha (default 0.05), otherwise "none".
# Hill estimates that iNEXT could not compute at the requested coverage are
# left out. A fit that fails, or a series with fewer than 3 years, gives NA.
#
# Output: output/indicators_v1/trends_<grid>_<level>.csv and
#         trend_direction_<grid>_<level>.csv

suppressPackageStartupMessages({ library(data.table); library(nlme) })
arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
GRID <- match.arg(arg("--grid", "1km"), c("100m", "1km"))
LEVEL <- as.numeric(arg("--level", "0.9"))
CI_FILE <- arg("--ci")
ALPHA <- as.numeric(arg("--alpha", "0.05"))
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

NO_FIT <- list(slope = NA_real_, p = NA_real_)
fit <- function(v, y, w = NULL) {
  if (length(v) < 3 || length(unique(v)) < 2) return(NO_FIT)
  co <- suppressWarnings(summary(if (is.null(w)) lm(v ~ y) else lm(v ~ y, weights = w)))$coefficients
  if (nrow(co) < 2) return(NO_FIT)
  list(slope = co[2, 1], p = co[2, 4])
}
# AR(1) errors on the year (gaps allowed); with se, variance proportional to se^2.
fit_gls <- function(v, y, se = NULL) {
  if (length(v) < 4 || length(unique(v)) < 2) return(NO_FIT)
  d <- data.frame(v = v, y = y, s2 = if (is.null(se)) 1 else se^2)
  g <- tryCatch(suppressWarnings(gls(v ~ y, data = d, correlation = corAR1(form = ~ y), method = "REML",
                                     weights = if (is.null(se)) NULL else varFixed(~ s2))),
                error = function(e) NULL)
  if (is.null(g)) return(NO_FIT)
  tt <- summary(g)$tTable
  list(slope = tt[2, 1], p = tt[2, 4])
}
direction <- function(slope, p) fifelse(is.na(p), NA_character_, fifelse(p <= ALPHA, fifelse(slope > 0, "up", "down"), "none"))

long <- merge(long, magv[, .(site_id, subset, stable_run_start, stable_run_end, has_usable_window)], by = c("site_id", "subset"), all.x = TRUE)
tr <- long[order(year), {
  f <- fit(value, year)
  st <- has_usable_window[1] %in% TRUE & !is.na(stable_run_start[1])
  fs <- if (st) { k <- year >= stable_run_start[1] & year <= stable_run_end[1]; fit(value[k], year[k]) } else NO_FIT
  okw <- !is.na(se)
  fw <- if (sum(okw) >= 3) fit(value[okw], year[okw], 1 / se[okw]^2) else NO_FIT
  fg <- fit_gls(value, year)
  fgw <- if (sum(okw) >= 4) fit_gls(value[okw], year[okw], se[okw]) else NO_FIT
  .(n_years = .N, first_year = min(year), last_year = max(year), mean_value = mean(value),
    slope = f$slope, p = f$p, slope_stable = fs$slope, p_stable = fs$p,
    slope_gls = fg$slope, p_gls = fg$p,
    n_years_with_se = sum(okw), slope_weighted = fw$slope, p_weighted = fw$p,
    slope_gls_weighted = fgw$slope, p_gls_weighted = fgw$p)
}, by = .(site_id, subset, indicator)]
tr[, `:=`(direction = direction(slope, p), direction_gls = direction(slope_gls, p_gls),
          direction_weighted = direction(slope_weighted, p_weighted),
          direction_gls_weighted = direction(slope_gls_weighted, p_gls_weighted), alpha = ALPHA)]
tag <- paste0(GRID, "_", LEVEL)
fwrite(tr, file.path(OUT, paste0("trends_", tag, ".csv")))
dir <- tr[!is.na(direction), .(combinations = .N, up = sum(direction == "up"), none = sum(direction == "none"),
                               down = sum(direction == "down")), by = .(subset, indicator)]
fwrite(dir, file.path(OUT, paste0("trend_direction_", tag, ".csv")))
cat(sprintf("Grid %s, coverage %.2f: %d trend rows written to %s\n", GRID, LEVEL, nrow(tr), file.path(OUT, paste0("trends_", tag, ".csv"))))
