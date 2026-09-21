#!/usr/bin/env Rscript
library(dplyr)

# Elbow threshold function
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

# 1. Load raw results (which now has corrected chao2 values)
raw_results_file <- "output/data_sufficiency/subgroup_analysis/raw_tests_results.csv"
if (!file.exists(raw_results_file)) {
  stop("raw_tests_results.csv not found!")
}
raw_res <- read.csv(raw_results_file, stringsAsFactors = FALSE)

# 2. Load the old precision audit to get mean_uncertainty_m and site_area_km2
old_precision_file <- "output/data_sufficiency/subgroup_analysis/all_tests_results_with_precision.csv"
if (!file.exists(old_precision_file)) {
  stop("all_tests_results_with_precision.csv not found!")
}
old_prec <- read.csv(old_precision_file, stringsAsFactors = FALSE) %>%
  select(continent, country, site_id, mean_uncertainty_m, site_area_km2) %>%
  distinct(continent, country, site_id, .keep_all = TRUE)

# 3. Join
results_full <- raw_res %>%
  left_join(old_prec, by = c("continent", "country", "site_id"))

# Calculate precision ratio and pass flags
DENSITY_THRESH <- 0.25
CHAO2_THRESH <- 0.70
SAC_SLOPE_THRESH <- 0.10
PRECISION_THRESH <- 0.10

# Calculate dynamic temporal threshold via elbow method
temporal_threshold <- elbow_threshold(results_full$temporal_distance)
cat(sprintf("Dynamic Temporal Threshold calculated via Elbow Method: %.4f\n", temporal_threshold))

results_full <- results_full %>%
  mutate(
    uncertainty_area_km2 = pi * ((mean_uncertainty_m / 1000)^2),
    precision_ratio = uncertainty_area_km2 / site_area_km2,
    
    pass_precision = !is.na(precision_ratio) & precision_ratio <= PRECISION_THRESH,
    pass_density = !is.na(density) & density >= DENSITY_THRESH,
    pass_chao2 = !is.na(chao2) & chao2 >= CHAO2_THRESH,
    pass_sac_slope = !is.na(sac_slope) & sac_slope <= SAC_SLOPE_THRESH,
    pass_temporal = !is.na(temporal_distance) & temporal_distance >= temporal_threshold,
    
    # FINAL PASS: 5 criteria conjoined (density, chao2, sac_slope, temporal, precision)
    pass_all = pass_precision & pass_density & pass_chao2 & pass_sac_slope & pass_temporal
  )

# Write output files
write.csv(results_full, "output/data_sufficiency/subgroup_analysis/all_tests_results_corrected.csv", row.names = FALSE)

# Sites passing corrected: only pass_all == TRUE, columns: continent, country, site_id, subset
sites_passing <- results_full %>%
  filter(pass_all) %>%
  select(continent, country, site_id, subset)

write.csv(sites_passing, "output/data_sufficiency/subgroup_analysis/sites_passing_corrected.csv", row.names = FALSE)

cat(sprintf("Saved corrected results. Passing site-subgroups: %d\n", nrow(sites_passing)))
