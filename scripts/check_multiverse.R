#!/usr/bin/env Rscript
# Robustness of the Hill-number trend verdicts across analysis choices (a
# "multiverse"): for each passing site x taxon combination, the verdict (up /
# none / down) under every combination of Hill order, grid, coverage level and
# trend model, and how often the verdicts agree.
#
# Usage (from the repository root, after summarise_indicators.R has been run
# for every grid and level used, with --ci for the 1 km grid):
#   Rscript scripts/check_multiverse.R                                  # 18 specifications
#   Rscript scripts/check_multiverse.R --models ols,wls,gls,gls_w --precision-rule median --alpha 0.05
#
# Specifications: q = 0, 1, 2 on each of
#   100 m and 1 km, at 90% and 95% coverage, for each model in --models;
#   weighted models (wls, gls_w) only where intervals exist (the 1 km grid).
# The default (--models ols,wls) gives 3 x (4 + 2) = 18 specifications.
# Verdicts are recomputed here from each slope and p value: "up"/"down" when
# p <= --alpha (default 0.05), otherwise "none".
# A combination is classed "always up", "always none", "always down",
# "up and down" (significant trends of both signs), "mixed" (a trend under some
# specifications and none under others) or "not all estimable" (some
# specification has no verdict, e.g. too few years reach the coverage level).
#
# Output: output/v1_checks/multiverse_matrix_<tag>.csv (one row per
# combination, one column per specification) and multiverse_per_combo_<tag>.csv,
# where <tag> is the precision rule, the alpha and any non-default models.

suppressPackageStartupMessages(library(data.table))
arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
MODELS <- strsplit(arg("--models", "ols,wls"), ",")[[1]]
stopifnot(all(MODELS %in% c("ols", "wls", "gls", "gls_w")))
RULE <- match.arg(arg("--precision-rule", "mean"), c("mean", "median", "d61"))
ALPHA <- as.numeric(arg("--alpha", "0.05"))
IN <- "output/indicators_v1"; OUT <- "output/v1_checks"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

COLS <- list(ols = c("slope", "p"), wls = c("slope_weighted", "p_weighted"),
             gls = c("slope_gls", "p_gls"), gls_w = c("slope_gls_weighted", "p_gls_weighted"))
verdict_of <- function(slope, p) fifelse(is.na(p), NA_character_, fifelse(p <= ALPHA, fifelse(slope > 0, "up", "down"), "none"))
combos <- unique(fread(file.path(IN, "ci_combinations.csv"))[precision_rule == RULE, .(site_id, subset)])

specs <- list()
for (g in c("100m", "1km")) for (lv in c("0.9", "0.95")) {
  f <- file.path(IN, sprintf("trends_%s_%s.csv", g, lv))
  if (!file.exists(f)) stop("Missing ", f, ": run summarise_indicators.R for this grid and level first")
  t <- merge(fread(f)[indicator %in% c("hill0", "hill1", "hill2")], combos, by = c("site_id", "subset"))
  for (m in MODELS) {
    if (m %in% c("wls", "gls_w") && g != "1km") next
    cc <- COLS[[m]]
    if (!all(cc %in% names(t))) stop(f, " has no ", cc[2], " column (re-run summarise_indicators.R, with --ci for weighted models)")
    specs[[length(specs) + 1]] <- t[, .(site_id, subset, spec = paste(indicator, g, lv, m),
                                        dir = verdict_of(get(cc[1]), get(cc[2])))]
  }
}
s <- rbindlist(specs)
n_spec <- uniqueN(s$spec)
wide <- dcast(s, site_id + subset ~ spec, value.var = "dir")
wide <- merge(combos, wide, by = c("site_id", "subset"), all.x = TRUE)

per <- melt(wide, id.vars = c("site_id", "subset"), variable.name = "spec", value.name = "dir")[, .(
  n_spec = n_spec, n_estimable = sum(!is.na(dir)),
  up = sum(dir %in% "up"), none = sum(dir %in% "none"), down = sum(dir %in% "down")), by = .(site_id, subset)]
per[, verdict := fifelse(n_estimable < n_spec, "not all estimable",
                 fifelse(up == n_spec, "always up", fifelse(down == n_spec, "always down",
                 fifelse(none == n_spec, "always none", fifelse(up > 0 & down > 0, "up and down", "mixed")))))]

tag <- paste0(RULE, "_a", ALPHA, if (!identical(MODELS, c("ols", "wls"))) paste0("_", paste(MODELS, collapse = "-")) else "")
fwrite(wide, file.path(OUT, sprintf("multiverse_matrix_%s.csv", tag)))
fwrite(per, file.path(OUT, sprintf("multiverse_per_combo_%s.csv", tag)))
cat(sprintf("Precision rule %s, models %s, alpha %.2f: %d combinations, %d specifications each.\n",
            RULE, paste(MODELS, collapse = ","), ALPHA, nrow(per), n_spec))
print(per[, .N, by = verdict][order(-N)])
