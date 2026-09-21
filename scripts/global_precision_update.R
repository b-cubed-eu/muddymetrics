# Global Precision Audit and Stats Update
library(dplyr)
library(sf)
library(parallel)
sf_use_s2(FALSE)

cat("Loading raw results...\n")
raw_results_file <- "output/data_sufficiency/subgroup_analysis/raw_tests_results.csv"
raw_results <- read.csv(raw_results_file)

# Helper functions
get_mean_uncertainty <- function(cont, country, site_id) {
  data_dir <- file.path("inst/extdata", paste0("ramsar_site_data_100m_", cont), country)
  cube_files <- list.files(data_dir, pattern = paste0("^", site_id, "(_|\\.|$)"), full.names = TRUE)
  if (length(cube_files) == 0) return(NA)
  tryCatch({
    header <- names(data.table::fread(cube_files[1], nrows = 0))
    uncertainty_cols <- c("coordinateuncertaintyinmeters", "mincoordinateuncertainty", "mincoordinateuncertaintyinmeters")
    existing_unc <- intersect(header, uncertainty_cols)
    if (length(existing_unc) > 0) {
      data <- data.table::fread(cube_files[1], select = existing_unc[1])
      return(mean(data[[existing_unc[1]]], na.rm = TRUE))
    }
    return(NA)
  }, error = function(e) return(NA))
}

get_site_area <- function(site_id) {
  wkt_files <- list.files("inst/extdata/ramsar_sites_wkt", pattern = paste0("^", site_id, "(_|\\.|$)"), recursive = TRUE, full.names = TRUE)
  if (length(wkt_files) == 0) return(NA)
  tryCatch({
    wkt_text <- paste(readLines(wkt_files[1], warn = FALSE), collapse = " ")
    poly <- st_make_valid(st_as_sfc(wkt_text, crs = 4326))
    return(as.numeric(st_area(poly)) / 1000000) # km^2
  }, error = function(e) return(NA))
}

# Unique sites for area calculation to avoid redundant WKT processing
unique_sites <- unique(raw_results$site_id)
cat(sprintf("Calculating areas for %d unique sites in parallel...\n", length(unique_sites)))

num_cores <- min(detectCores() - 1, 10)
cl_area <- makePSOCKcluster(num_cores)
clusterExport(cl_area, c("get_site_area"))
clusterEvalQ(cl_area, { library(sf); sf_use_s2(FALSE) })
site_areas_list <- parLapply(cl_area, unique_sites, get_site_area)
stopCluster(cl_area)

site_areas <- setNames(site_areas_list, unique_sites)

# Parallel uncertainty extraction
unique_sites_df <- raw_results %>%
  select(continent, country, site_id) %>%
  distinct()

cat(sprintf("Extracting uncertainty for %d unique sites in parallel...\n", nrow(unique_sites_df)))
cl <- makePSOCKcluster(num_cores)
clusterExport(cl, c("get_mean_uncertainty", "unique_sites_df"), envir = environment())
clusterEvalQ(cl, { library(data.table) })

uncertainties_list <- parLapply(cl, seq_len(nrow(unique_sites_df)), function(i) {
  get_mean_uncertainty(unique_sites_df$continent[i], unique_sites_df$country[i], unique_sites_df$site_id[i])
})
stopCluster(cl)

unique_sites_df$mean_uncertainty_m <- unlist(uncertainties_list)

# Combine and apply 0.1 cutoff
cat("Combining results and applying 0.1 cutoff...\n")
# Thresholds from compute_thresholds.R
DENSITY_THRESH <- 0.25
CHAO2_THRESH <- 0.7
SAC_SLOPE_THRESH <- 0.1
TEMPORAL_THRESH <- 0.25
PRECISION_THRESH <- 0.1

results_full <- raw_results %>%
  left_join(unique_sites_df, by = c("continent", "country", "site_id")) %>%
  mutate(
    site_area_km2 = as.numeric(site_areas[site_id]),
    uncertainty_area_km2 = pi * ((mean_uncertainty_m / 1000)^2),
    precision_ratio = uncertainty_area_km2 / site_area_km2,
    # New spatial pass logic: Precision Ratio <= 0.1
    pass_precision = !is.na(precision_ratio) & precision_ratio <= PRECISION_THRESH,
    # Standard passes (handling NAs)
    pass_density = !is.na(density) & density >= DENSITY_THRESH,
    pass_chao2 = !is.na(chao2) & chao2 >= CHAO2_THRESH,
    pass_sac_slope = !is.na(sac_slope) & sac_slope <= SAC_SLOPE_THRESH,
    pass_temporal = !is.na(temporal_distance) & temporal_distance >= TEMPORAL_THRESH,
    # FINAL PASS
    pass_all = pass_density & pass_chao2 & pass_sac_slope & pass_temporal & pass_precision
  )

# Update Stats
update_stats <- function(df, continent_name) {
  df %>%
    summarise(
      pass_density_count = sum(pass_density, na.rm = TRUE),
      pass_density_pct = mean(pass_density, na.rm = TRUE) * 100,
      pass_chao2_count = sum(pass_chao2, na.rm = TRUE),
      pass_chao2_pct = mean(pass_chao2, na.rm = TRUE) * 100,
      pass_sac_slope_count = sum(pass_sac_slope, na.rm = TRUE),
      pass_sac_slope_pct = mean(pass_sac_slope, na.rm = TRUE) * 100,
      pass_temporal_count = sum(pass_temporal, na.rm = TRUE),
      pass_temporal_pct = mean(pass_temporal, na.rm = TRUE) * 100,
      pass_precision_count = sum(pass_precision, na.rm = TRUE),
      pass_precision_pct = mean(pass_precision, na.rm = TRUE) * 100,
      pass_all_count = sum(pass_all, na.rm = TRUE),
      pass_all_pct = mean(pass_all, na.rm = TRUE) * 100,
      continent = continent_name,
      total_site_subsets = n()
    )
}

overall_summary <- update_stats(results_full, "All")
cont_summary <- results_full %>% group_by(continent) %>% group_map(~update_stats(.x, .y$continent[[1]])) %>% bind_rows()
final_stats <- bind_rows(overall_summary, cont_summary)

write.csv(final_stats, "output/data_sufficiency/subgroup_analysis/summary_stats_precision.csv", row.names = FALSE)
write.csv(results_full %>% filter(pass_all), "output/data_sufficiency/subgroup_analysis/sites_passing_all_precision.csv", row.names = FALSE)
write.csv(results_full, "output/data_sufficiency/subgroup_analysis/all_tests_results_with_precision.csv", row.names = FALSE)

cat("✓ New statistics and passing list saved.\n")
print(overall_summary)
