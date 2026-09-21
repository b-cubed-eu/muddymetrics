#!/usr/bin/env Rscript

library(dplyr)

DENSITY_THRESHOLD <- 0.25
CHAO2_THRESHOLD <- 0.7
SAC_SLOPE_THRESHOLD <- 0.1

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

out_dir <- "output/data_sufficiency/subgroup_analysis"
in_file <- file.path(out_dir, "raw_tests_results.csv")

if (!file.exists(in_file)) {
  cat("Raw results file not found.\n")
  quit(status=1)
}

res_df <- read.csv(in_file, stringsAsFactors = FALSE)

# Compute elbow thresholds
spatial_threshold <- elbow_threshold(res_df$spatial_distance)
temporal_threshold <- elbow_threshold(res_df$temporal_distance)

cat(sprintf("\nDetermined thresholds via Elbow Method:\nSpatial Threshold: %.4f\nTemporal Threshold: %.4f\n\n", spatial_threshold, temporal_threshold))

res_df$pass_density <- as.integer(!is.na(res_df$density) & res_df$density >= DENSITY_THRESHOLD)
res_df$pass_chao2 <- as.integer(!is.na(res_df$chao2) & res_df$chao2 >= CHAO2_THRESHOLD)
res_df$pass_sac_slope <- as.integer(!is.na(res_df$sac_slope) & res_df$sac_slope <= SAC_SLOPE_THRESHOLD)
res_df$pass_temporal <- as.integer(!is.na(res_df$temporal_distance) & res_df$temporal_distance >= temporal_threshold)
res_df$pass_spatial <- as.integer(!is.na(res_df$spatial_distance) & res_df$spatial_distance >= spatial_threshold)

res_df$pass_all <- as.integer(res_df$pass_density == 1 & res_df$pass_chao2 == 1 & res_df$pass_sac_slope == 1 & res_df$pass_temporal == 1 & res_df$pass_spatial == 1)
res_df$pass_all_except_spatial <- as.integer(res_df$pass_density == 1 & res_df$pass_chao2 == 1 & res_df$pass_sac_slope == 1 & res_df$pass_temporal == 1)
res_df$pass_all_except_temporal <- as.integer(res_df$pass_density == 1 & res_df$pass_chao2 == 1 & res_df$pass_sac_slope == 1 & res_df$pass_spatial == 1)

write.csv(res_df, file.path(out_dir, "all_tests_results.csv"), row.names = FALSE)

write.csv(res_df[res_df$pass_all == 1, c("site_id", "continent", "country", "subset")], 
          file.path(out_dir, "sites_passing_all.csv"), row.names = FALSE)
write.csv(res_df[res_df$pass_all_except_spatial == 1, c("site_id", "continent", "country", "subset")], 
          file.path(out_dir, "sites_passing_all_except_spatial.csv"), row.names = FALSE)
write.csv(res_df[res_df$pass_all_except_temporal == 1, c("site_id", "continent", "country", "subset")], 
          file.path(out_dir, "sites_passing_all_except_temporal.csv"), row.names = FALSE)
          
test_cols <- c("pass_density", "pass_chao2", "pass_sac_slope", "pass_temporal", "pass_spatial", "pass_all")

overall_summary <- res_df %>%
  summarise(across(all_of(test_cols), list(
    count = ~sum(., na.rm = TRUE),
    pct = ~sum(., na.rm = TRUE) / n() * 100
  ))) %>%
  mutate(continent = "All", total_site_subsets = nrow(res_df))
  
cont_summary <- res_df %>%
  group_by(continent) %>%
  summarise(
    total_site_subsets = n(),
    across(all_of(test_cols), list(
      count = ~sum(., na.rm = TRUE),
      pct = ~sum(., na.rm = TRUE) / n() * 100
    )),
    .groups = "drop"
  )
  
summary_stats <- bind_rows(overall_summary, cont_summary)
write.csv(summary_stats, file.path(out_dir, "summary_stats.csv"), row.names = FALSE)

cat("Threshold and stats calculation complete. Check output directory.\n")
