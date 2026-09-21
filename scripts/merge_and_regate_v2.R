#!/usr/bin/env Rscript
# =============================================================================
# merge_and_regate_v2.R
#
# Joins the recomputed temporal distances and occurrence totals onto the
# published results, fixes the density units, re-derives the dynamic elbow
# threshold, and re-applies the five gates.
#
# NON-DESTRUCTIVE: reads all_tests_results_corrected.csv and
# temporal_only_results.csv; writes only *_v2 files.
#
# Columns carried over UNCHANGED from the published results:
#   chao2, sac_slope, mean_uncertainty_m, site_area_km2,
#   uncertainty_area_km2, precision_ratio
# Columns recomputed:
#   temporal_distance (cumsum removed), density (units fixed)
#
# GATE 0 IS A HARD REGRESSION TEST. If the reproduced cumsum statistic does not
# match the published temporal_distance column, the script stops. That check is
# what licenses trusting the corrected value.
# =============================================================================
library(dplyr)

OUT_DIR <- "output/data_sufficiency/subgroup_analysis"

CHAO2_THRESH     <- 0.70
SAC_SLOPE_THRESH <- 0.10
PRECISION_THRESH <- 0.10
DENSITY_THRESH   <- 0.25   # occurrences per km^2, Troia & McManamay (2016)

# Identical to global_precision_update_corrected.R and compute_thresholds.R.
elbow_threshold <- function(distances) {
  dist_clean <- sort(distances[!is.na(distances)])
  if (length(dist_clean) < 10) return(0.25)
  n <- length(dist_clean)
  x <- (1:n) / n
  y <- (dist_clean - min(dist_clean)) / (max(dist_clean) - min(dist_clean))
  dist_to_line <- abs(x - y) / sqrt(2)
  elbow_idx <- which.max(dist_to_line)
  if (length(elbow_idx) == 0) return(0.25)
  return(dist_clean[elbow_idx])
}

orig <- read.csv(file.path(OUT_DIR, "all_tests_results_corrected.csv"), stringsAsFactors = FALSE)
new  <- read.csv(file.path(OUT_DIR, "temporal_only_results.csv"), stringsAsFactors = FALSE)

cat(sprintf("published rows: %d | recomputed rows: %d\n", nrow(orig), nrow(new)))

j <- orig %>%
  inner_join(new %>% select(site_id, subset, temporal_distance_cumsum,
                            temporal_distance_annual, total_occurrences),
             by = c("site_id", "subset"))
cat(sprintf("joined rows: %d\n", nrow(j)))

# ---------------------------------------------------------------------------
# GATE 0 — regression test
# ---------------------------------------------------------------------------
cat("\n=============== GATE 0: REGRESSION TEST ===============\n")
both <- j %>% filter(!is.na(temporal_distance) & !is.na(temporal_distance_cumsum))
agree <- sum(abs(both$temporal_distance - both$temporal_distance_cumsum) < 1e-9)
pct <- 100 * agree / nrow(both)
cat(sprintf("reproduced cumsum statistic matches published: %d / %d (%.2f%%)\n", agree, nrow(both), pct))

if (pct < 99) {
  cat("\nFAIL. The worker is not reproducing the original computation.\n")
  cat("Do not trust temporal_distance_annual from this run. Investigate before proceeding.\n")
  mism <- both %>% filter(abs(temporal_distance - temporal_distance_cumsum) >= 1e-9) %>%
    select(site_id, subset, published = temporal_distance, reproduced = temporal_distance_cumsum)
  print(utils::head(as.data.frame(mism), 20))
  write.csv(mism, file.path(OUT_DIR, "regression_mismatches_v2.csv"), row.names = FALSE)
  stop("Regression test failed.")
}
cat("PASS. The worker reproduces the original computation; corrected values are trustworthy.\n")

# ---------------------------------------------------------------------------
# Elbow thresholds
# ---------------------------------------------------------------------------
cat("\n=============== ELBOW THRESHOLDS ===============\n")
elbow_old <- elbow_threshold(orig$temporal_distance)
elbow_new <- elbow_threshold(j$temporal_distance_annual)
cat(sprintf("elbow on published (broken) distances : %.4f   [must be 0.2500]\n", elbow_old))
cat(sprintf("elbow on corrected distances          : %.4f\n", elbow_new))
old_pass_repro <- sum(!is.na(orig$temporal_distance) & orig$temporal_distance >= elbow_old)
cat(sprintf("published pass_temporal reproduced    : %d   [must be 1357]\n", old_pass_repro))
if (abs(elbow_old - 0.25) > 1e-9 || old_pass_repro != 1357) {
  stop("Elbow function does not reproduce the published gate. Stop.")
}
cat("PASS.\n")

# ---------------------------------------------------------------------------
# Density fix
# ---------------------------------------------------------------------------
cat("\n=============== DENSITY UNIT FIX ===============\n")
j <- j %>%
  mutate(
    density_published = density,
    density_km2 = ifelse(!is.na(site_area_km2) & site_area_km2 > 0,
                         total_occurrences / site_area_km2, NA_real_),
    density_ratio = density_published / density_km2
  )
dr <- j$density_ratio[is.finite(j$density_ratio)]
cat(sprintf("median published/corrected density ratio: %.0f\n", median(dr, na.rm = TRUE)))
cat("Expected 4,000-12,500 — one square degree is 12,321 km^2 at the equator,\n")
cat("scaling with cos(latitude). A ratio near 1 would mean the diagnosis was wrong.\n")
if (median(dr, na.rm = TRUE) < 3000) {
  cat("\nWARNING: ratio is not consistent with a square-degree/km^2 confusion.\n")
  cat("Re-examine before using the corrected density.\n")
}
cat(sprintf("corrected density: median %.2f, 5th pct %.4f occurrences/km2\n",
            median(j$density_km2, na.rm = TRUE), quantile(j$density_km2, 0.05, na.rm = TRUE)))

# ---------------------------------------------------------------------------
# Re-apply gates
# ---------------------------------------------------------------------------
res <- j %>%
  mutate(
    temporal_distance_v2 = temporal_distance_annual,
    pass_precision = !is.na(precision_ratio) & precision_ratio <= PRECISION_THRESH,
    pass_chao2     = !is.na(chao2) & chao2 >= CHAO2_THRESH,
    pass_sac_slope = !is.na(sac_slope) & sac_slope <= SAC_SLOPE_THRESH,
    pass_temporal  = !is.na(temporal_distance_v2) & temporal_distance_v2 >= elbow_new,
    pass_density_published = !is.na(density_published) & density_published >= DENSITY_THRESH,
    pass_density   = !is.na(density_km2) & density_km2 >= DENSITY_THRESH,
    # headline: all five, with both fixes applied
    pass_all = pass_precision & pass_density & pass_chao2 & pass_sac_slope & pass_temporal,
    # isolates the temporal fix alone, density left as published
    pass_all_temporal_fix_only = pass_precision & pass_density_published &
                                 pass_chao2 & pass_sac_slope & pass_temporal
  )

cat("\n=============== CASCADE ===============\n")
cat(sprintf("total combinations              : %d\n", nrow(res)))
cat(sprintf("pass density (published units)  : %d\n", sum(res$pass_density_published)))
cat(sprintf("pass density (corrected km2)    : %d\n", sum(res$pass_density)))
cat(sprintf("pass chao2                      : %d\n", sum(res$pass_chao2)))
cat(sprintf("pass sac slope                  : %d\n", sum(res$pass_sac_slope)))
cat(sprintf("pass precision                  : %d\n", sum(res$pass_precision)))
cat(sprintf("pass temporal (corrected)       : %d\n", sum(res$pass_temporal)))
cat(sprintf("\nPASS ALL, temporal fix only     : %d combos\n", sum(res$pass_all_temporal_fix_only)))
cat(sprintf("PASS ALL, both fixes            : %d combos, %d unique sites\n",
            sum(res$pass_all), length(unique(res$site_id[res$pass_all]))))
cat(sprintf("(published was 77 combos)\n"))

cat("\n=============== TAXONOMIC COMPOSITION ===============\n")
old_pass <- orig %>% filter(as.logical(pass_all))
comp <- full_join(old_pass %>% count(subset, name = "published"),
                  res %>% filter(pass_all) %>% count(subset, name = "corrected"),
                  by = "subset") %>%
  mutate(across(c(published, corrected), ~ifelse(is.na(.), 0L, .)),
         published_pct = round(100 * published / sum(published), 1),
         corrected_pct = round(100 * corrected / sum(corrected), 1))
print(as.data.frame(comp))

cat("\n=============== PART 3 VIGNETTE SITES ===============\n")
print(as.data.frame(res %>%
  filter(site_id %in% c("site_786", "site_67", "site_879"), subset %in% c("full", "aves")) %>%
  select(site_id, subset, temporal_distance_v2, density_km2, pass_temporal, pass_density, pass_all)))

cat("\n=============== UNCHANGED-COLUMN CHECK ===============\n")
# chao2, sac_slope and precision are carried through the join, never recomputed,
# so they must be identical by construction. This verifies the join did not
# reorder or duplicate anything.
chk <- res %>%
  select(site_id, subset, chao2, sac_slope, precision_ratio) %>%
  inner_join(orig %>% select(site_id, subset,
                             chao2_o = chao2, sac_o = sac_slope, prec_o = precision_ratio),
             by = c("site_id", "subset"))
cat(sprintf("rows after re-join   : %d (must equal %d)\n", nrow(chk), nrow(res)))
cat(sprintf("chao2 identical      : %d / %d\n",
            sum(chk$chao2 == chk$chao2_o | (is.na(chk$chao2) & is.na(chk$chao2_o))), nrow(chk)))
cat(sprintf("sac_slope identical  : %d / %d\n",
            sum(chk$sac_slope == chk$sac_o | (is.na(chk$sac_slope) & is.na(chk$sac_o))), nrow(chk)))
cat(sprintf("precision identical  : %d / %d\n",
            sum(chk$precision_ratio == chk$prec_o | (is.na(chk$precision_ratio) & is.na(chk$prec_o))), nrow(chk)))
if (nrow(chk) != nrow(res)) stop("Join integrity check failed - row count changed.")

lost   <- res %>% filter(!pass_all) %>% semi_join(old_pass, by = c("site_id", "subset"))
gained <- res %>% filter(pass_all) %>% anti_join(old_pass, by = c("site_id", "subset"))
cat(sprintf("\nlost   (passed before, fail now): %d\n", nrow(lost)))
cat(sprintf("gained (failed before, pass now): %d\n", nrow(gained)))

write.csv(res,    file.path(OUT_DIR, "all_tests_results_v2.csv"), row.names = FALSE)
write.csv(res %>% filter(pass_all), file.path(OUT_DIR, "sites_passing_v2.csv"), row.names = FALSE)
write.csv(lost,   file.path(OUT_DIR, "temporal_fix_lost_v2.csv"), row.names = FALSE)
write.csv(gained, file.path(OUT_DIR, "temporal_fix_gained_v2.csv"), row.names = FALSE)

cat("\nWrote all_tests_results_v2.csv, sites_passing_v2.csv,\n")
cat("      temporal_fix_lost_v2.csv, temporal_fix_gained_v2.csv\n")
