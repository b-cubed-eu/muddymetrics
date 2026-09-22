#!/usr/bin/env Rscript
# After a full run of run_sufficiency.R, compare the new raw values with the
# D6.1 published results and the August 2026 (v2) correction.
#
# Usage (from the repository root):
#   Rscript scripts/check_sufficiency_against_d61.R
#
# Reads:
#   output/sufficiency_v1/sufficiency_raw.csv
#   output/data_sufficiency/subgroup_analysis/all_tests_results_v2.csv
#
# What to expect (full run, 6,301 non-herps combinations):
#   EXACT  site area, total occurrences, density, D6.1 precision definition,
#          Chao2 by key (except where D6.1 kept an old value because Chao2 is
#          undefined), annual coverage (b3gbi rule)
#   DIFFER binary temporal statistic: v2 used b3gbi 0.9.x, whose annual
#          richness is wrong where species keys are alphanumeric
#   DIFFER Chao2 by name: only where a species name has more than one key
# Any EXACT check that fails at full scale should be investigated before the
# results are used.

suppressPackageStartupMessages(library(dplyr))
new <- read.csv("output/sufficiency_v1/sufficiency_raw.csv", stringsAsFactors = FALSE)
v2  <- read.csv("output/data_sufficiency/subgroup_analysis/all_tests_results_v2.csv", stringsAsFactors = FALSE)

j <- new %>% filter(status == "computed", subset != "herps") %>%
  inner_join(v2 %>% select(site_id, subset, v2_area = site_area_km2, v2_density = density_km2,
                           v2_total_occ = total_occurrences, v2_td_annual = temporal_distance_annual,
                           pub_chao2 = chao2, pub_sac = sac_slope, pub_prec = precision_ratio),
             by = c("site_id", "subset"))

same <- function(a, b, tol = 1e-9) ifelse(is.na(a) & is.na(b), TRUE,
                                   ifelse(is.na(a) | is.na(b), FALSE, abs(a - b) <= tol * pmax(1, abs(b))))
check <- function(label, ok, expect) {
  status <- if (all(ok)) "EXACT" else sprintf("differs in %d", sum(!ok))
  flag <- if (expect == "exact" && !all(ok)) "  <-- INVESTIGATE" else ""
  cat(sprintf("  %-50s %6d / %-6d %s%s\n", label, sum(ok), length(ok), status, flag))
}

cat(sprintf("\nNew run: %d sites, %d site x taxon rows, %d computed\n",
            n_distinct(new$site_id), nrow(new), sum(new$status == "computed")))
cat(sprintf("Non-herps combinations present in both new and v2: %d\n\n", nrow(j)))
cat("Expected EXACT:\n")
check("site area (km2)", same(j$site_area_km2, j$v2_area), "exact")
check("total occurrences", same(j$n_occurrences, j$v2_total_occ), "exact")
# Site areas of very small sites can differ in the 9th significant digit between
# sf/GEOS/lwgeom builds, so derived values are compared to 1 part in a million.
check("density (occ/km2)", same(j$density_km2, j$v2_density, 1e-6), "exact")
check("precision, D6.1 definition", same(j$precision_ratio_d61, j$pub_prec, 1e-6), "exact")
fb <- is.na(j$chao2_completeness_bykey) & !is.na(j$pub_chao2)   # the old fallback, see below
check("Chao2 completeness by key (excluding old fallback)", same(j$chao2_completeness_bykey, j$pub_chao2)[!fb], "exact")
check("annual coverage deficit (b3gbi rule)", same(j$coverage_deficit_b3gbi, j$pub_sac, 1e-6), "exact")
cat("\nExpected to differ (explained):\n")
check("binary temporal statistic vs v2", same(j$temporal_binary, j$v2_td_annual), "differ")
ch <- same(j$chao2_completeness, j$pub_chao2)
check("Chao2 completeness by name", ch, "differ")
cat(sprintf("    differences at combinations with split species names: %d of %d\n",
            sum(!ch & j$n_names_multi_key > 0), sum(!ch)))
cat(sprintf("    published Chao2 present where it is undefined (the old fallback): %d\n", sum(fb)))

missing_in_new <- v2 %>% filter(site_id %in% new$site_id) %>% select(site_id, subset) %>%
  anti_join(new %>% filter(status == "computed"), by = c("site_id", "subset"))
cat(sprintf("\nv2 combinations (at sites in the new run) not computed in the new run: %d (expected 0)\n", nrow(missing_in_new)))
if (file.exists("output/sufficiency_v1/failures.csv"))
  cat("NOTE: output/sufficiency_v1/failures.csv exists - some sites failed; see output/sufficiency_v1/logs/\n")
