#!/usr/bin/env Rscript
# Calculate density by class for multiple thresholds - simplified version
# Uses pre-calculated density_by_class_1990.csv and applies different thresholds

library(dplyr)

# Read the pre-calculated data (already has 0.25 threshold applied)
density <- read.csv('output/data_sufficiency/density_by_class_1990.csv')

# Create directories
threshold_names <- c("liberal_0.10", "standard_0.25", "conservative_0.50")
thresholds <- c(0.10, 0.25, 0.50)

for(tn in threshold_names){
  dir.create(paste0("output/data_sufficiency/", tn), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("output/data_sufficiency/", tn, "/by_class"), showWarnings = FALSE)
  dir.create(paste0("output/data_sufficiency/", tn, "/by_continent"), showWarnings = FALSE)
}

# Process each threshold
for(idx in 1:length(thresholds)){
  thresh <- thresholds[idx]
  thresh_name <- threshold_names[idx]
  
  cat("Processing threshold:", thresh, "\n")
  
  # Recalculate passes for this threshold
  results <- density
  results$passes <- as.numeric(results$density >= thresh)
  results$threshold_value <- thresh
  results$threshold_name <- thresh_name
  
  # Save full results
  write.csv(results, paste0('output/data_sufficiency/', thresh_name, '/density_by_class.csv'), row.names=FALSE)
  
  # Summary by class and continent
  summary <- results %>%
    group_by(continent, class) %>%
    summarise(
      total = n(),
      passing = sum(passes, na.rm=T),
      pass_rate = round(passing/total*100, 1),
      .groups = 'drop'
    )
  write.csv(summary, paste0('output/data_sufficiency/', thresh_name, '/by_class/summary.csv'), row.names=FALSE)
  
  # Summary by continent
  cont_summary <- results %>%
    group_by(continent) %>%
    summarise(
      total_sites = n_distinct(site_id),
      total_combos = n(),
      passing_combos = sum(passes, na.rm=T),
      pass_rate = round(passing_combos/total_combos*100, 1),
      .groups = 'drop'
    )
  write.csv(cont_summary, paste0('output/data_sufficiency/', thresh_name, '/by_continent/summary.csv'), row.names=FALSE)
  
  # Sites with >=1 passing class
  sites_pass <- results %>%
    group_by(continent, site_id, country, area_km2) %>%
    summarise(passes_any = as.numeric(any(passes == 1)), .groups = 'drop') %>%
    filter(passes_any == 1)
  write.csv(sites_pass, paste0('output/data_sufficiency/', thresh_name, '/sites_passing.csv'), row.names=FALSE)
  
  cat("  Sites passing:", nrow(sites_pass), "\n")
}

cat("\nDone!\n")
