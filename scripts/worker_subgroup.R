#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: worker_subgroup.R <continent> <country> <site_id> <out_file>\n")
  quit(status = 1)
}

cont <- args[1]
country <- args[2]
site_id <- args[3]
out_file <- args[4]

library(b3gbi)
library(dplyr)
library(sf)

DENSITY_THRESHOLD <- 0.25
CHAO2_THRESHOLD <- 0.7
SAC_SLOPE_THRESHOLD <- 0.1

SUBSETS <- c("full", "aves", "mammalia", "herps", "insects", "plants")

get_taxa_filter <- function(subset) {
  switch(subset,
    "full" = NULL,
    "aves" = function(df) df$class == "Aves",
    "mammalia" = function(df) df$class == "Mammalia",
    "herps" = function(df) df$class %in% c("Amphibia", "Reptilia"),
    "insects" = function(df) df$class == "Insecta",
    "plants" = function(df) df$kingdom == "Plantae",
    stop("Unknown subset: ", subset)
  )
}

# =============================================================================
# SUPERSEDED - DO NOT USE. Retained only for provenance.
# -----------------------------------------------------------------------------
# This is the ORIGINAL, INCORRECT Chao2 implementation from the first pipeline
# run. It did NOT produce any number reported in the deliverable.
#
# Two errors:
#   1. Wrong sampling unit. `table(valid_species)` counts CUBE ROWS per species,
#      i.e. species x cell x year combinations. The correct unit is the grid
#      cell, with years pooled, built as a species x cell incidence matrix.
#   2. No bias correction. The estimator below omits the (T-1)/T factor of the
#      bias-corrected Chao2.
#
# These are not cosmetic. Checked against 14 sites near the 0.70 sufficiency
# gate, this version differs from the correct one by up to 0.023 and would flip
# 4 of those 14 sites across the pass/fail threshold.
#
# THE CORRECT IMPLEMENTATION, and the one that produced every reported value,
# is `calculate_chao2()` in R/chao2_sac_functions.R, which builds a proper
# species x grid-cell incidence matrix and uses vegan::specpool(). Verified:
# 10 of 10 randomly sampled sites reproduce the reported figures exactly.
#
# If this worker is ever re-run, replace the call at the "Main Logic" section
# below with the R/chao2_sac_functions.R version before trusting the output.
# =============================================================================
calculate_chao2_DEPRECATED_ROWCOUNT <- function(data) {
  if (is.null(data) || nrow(data) < 10) return(NA_real_)
  valid_species <- data$specieskey[!is.na(data$specieskey)]
  if (length(valid_species) == 0) return(NA_real_)
  
  sample_counts <- table(valid_species)
  f1 <- sum(sample_counts == 1)
  f2 <- sum(sample_counts == 2)
  S_obs <- length(unique(valid_species))
  
  if (f2 > 0) {
    chao2_est <- S_obs + (f1^2) / (2 * f2)
  } else if (f1 > 0) {
    chao2_est <- S_obs + (f1 * (f1 - 1)) / 2
  } else {
    chao2_est <- S_obs
  }
  
  if (chao2_est > 0) return(S_obs / chao2_est)
  return(NA_real_)
}

calculate_jaccard <- function(occ_vec, rich_vec) {
  if (length(occ_vec) < 3) return(NA_real_)
  
  occ_norm <- (occ_vec - min(occ_vec)) / (max(occ_vec) - min(occ_vec) + 1e-10)
  rich_norm <- (rich_vec - min(rich_vec)) / (max(rich_vec) - min(rich_vec) + 1e-10)
  
  occ_threshold <- mean(occ_norm)
  rich_threshold <- mean(rich_norm)
  
  occ_set <- which(occ_norm > occ_threshold)
  rich_set <- which(rich_norm > rich_threshold)
  
  intersection <- length(intersect(occ_set, rich_set))
  union_set <- length(union(occ_set, rich_set))
  
  if (union_set == 0) return(NA_real_)
  jaccard_sim <- intersection / union_set
  return(1 - jaccard_sim)
}

get_site_area <- function(country, site_id) {
  wkt_dir <- file.path("inst/extdata", "ramsar_sites_wkt", country)
  if (!dir.exists(wkt_dir)) return(NA_real_)
  
  wkt_files <- list.files(wkt_dir, pattern = "\\.wkt$", full.names = FALSE)
  site_num <- strsplit(site_id, "_")[[1]][2] 
  
  wkt_file <- NULL
  for (wf in wkt_files) {
    if (grepl(paste0("^", site_num, "(_|\\.)"), wf) || grepl(paste0("^site_", site_num, "(_|\\.)"), wf) || wf == paste0(site_num, ".wkt")) {
      wkt_file <- file.path(wkt_dir, wf)
      break
    }
  }
  
  if (is.null(wkt_file) || !file.exists(wkt_file)) return(NA_real_)
  
  area <- tryCatch({
    wkt <- paste(readLines(wkt_file, warn = FALSE), collapse = "\n")
    poly <- sf::st_as_sfc(wkt)
    as.numeric(sf::st_area(poly) |> units::set_units("km^2"))
  }, error = function(e) NA_real_)
  return(area)
}

# Main Logic
processed_records <- tryCatch(read.csv(out_file, stringsAsFactors = FALSE), error = function(e) data.frame())
processed_ids <- if(nrow(processed_records) > 0) paste(processed_records$site_id, processed_records$subset, sep="_") else c()

data_dir <- file.path("inst/extdata", paste0("ramsar_site_data_100m_", cont))
cube_dir <- file.path(data_dir, country)
cube_files <- list.files(cube_dir, pattern = "\\.csv$", full.names = TRUE)

cube_path <- NULL
for (f in cube_files) {
  if (grepl(site_id, basename(f))) {
    cube_path <- f
    break
  }
}

if (is.null(cube_path)) quit(status = 0)

area <- get_site_area(country, site_id)

data <- tryCatch(read.csv(cube_path, stringsAsFactors = FALSE), error = function(e) NULL)
if (is.null(data) || nrow(data) < 10) quit(status = 0)
data <- data[data$year >= 1990, ]
if (nrow(data) < 10) quit(status = 0)

# temp_cube_file <- paste0("temp_worker_", site_id, ".csv") # No longer needed

for (subset in SUBSETS) {
  if (paste(site_id, subset, sep="_") %in% processed_ids) next
  
  taxa_filter <- get_taxa_filter(subset)
  if (is.null(taxa_filter)) {
    subset_data <- data
  } else {
    subset_data <- tryCatch(data[taxa_filter(data), ], error = function(e) data[FALSE, ])
  }
  
  if (nrow(subset_data) < 10) next
  
  density_val <- if (!is.na(area) && area > 0) sum(subset_data$occurrences, na.rm = TRUE) / area else NA_real_
  # SUPERSEDED CALL - see the deprecation notice on the Chao2 function above.
  # The original call here used the incorrect row-count implementation. It is
  # deliberately disabled rather than silently repointed, because re-running this
  # worker would overwrite raw_tests_results.csv, from which the reported
  # all_tests_results_corrected.csv is derived.
  #
  # To re-enable, source R/chao2_sac_functions.R and call its calculate_chao2(),
  # which takes a cube (or path) and returns a list - use $completeness here,
  # not the bare return value:
  #   chao2_val <- calculate_chao2(cube_path, shapefilepath)$completeness
  stop("Chao2 in this worker is superseded. See the notice above and use ",
       "R/chao2_sac_functions.R::calculate_chao2() before re-running.")
  chao2_val <- NA_real_
  
  cube <- tryCatch({
    b3gbi::process_cube(subset_data, cols_cellCode = "mgrscellcode", cols_speciesKey = "specieskey", cols_scientificName = "species")
  }, error = function(e) NULL)
  
  sac_slope <- NA_real_
  temporal_dist <- NA_real_
  spatial_dist <- NA_real_
  
  if (!is.null(cube)) {
    comp_ts <- tryCatch(as.data.frame(b3gbi::completeness_ts(cube)[["data"]]), error = function(e) NULL)
    if (!is.null(comp_ts) && "diversity_val" %in% names(comp_ts)) {
      mean_comp <- mean(comp_ts$diversity_val, na.rm = TRUE)
      if (!is.na(mean_comp)) sac_slope <- 1 - mean_comp
    }
    
    tot_occ_ts <- tryCatch(as.data.frame(b3gbi::total_occ_ts(cube)[["data"]]), error = function(e) NULL)
    obs_rich_ts <- tryCatch(as.data.frame(b3gbi::obs_richness_ts(cube)[["data"]]), error = function(e) NULL)
    if (!is.null(tot_occ_ts) && !is.null(obs_rich_ts) && nrow(tot_occ_ts) >= 3) {
      merged_ts <- merge(tot_occ_ts, obs_rich_ts, by = "year")
      merged_ts <- merged_ts[order(merged_ts$year), ]
      if (nrow(merged_ts) >= 3) {
        occ_cum <- cumsum(merged_ts$diversity_val.x)
        rich_cum <- cumsum(merged_ts$diversity_val.y)
        temporal_dist <- calculate_jaccard(occ_cum, rich_cum)
      }
    }
    
    tot_occ_map <- tryCatch(as.data.frame(b3gbi::total_occ_map(cube, cell_size="auto")[["data"]]), error = function(e) NULL)
    obs_rich_map <- tryCatch(as.data.frame(b3gbi::obs_richness_map(cube, cell_size="auto")[["data"]]), error = function(e) NULL)
    if (!is.null(tot_occ_map) && !is.null(obs_rich_map) && nrow(tot_occ_map) >= 3) {
      merged_map <- merge(tot_occ_map, obs_rich_map, by = "cellCode")
      if (nrow(merged_map) >= 3) {
        spatial_dist <- calculate_jaccard(merged_map$diversity_val.x, merged_map$diversity_val.y)
      }
    }
  }
  
  # Note: writing without locking, but master ensures synchronous workers
  cat(sprintf("%s,\"%s\",%s,%s,%s,%s,%s,%s,%s\n", 
              cont, country, site_id, subset, 
              ifelse(is.na(density_val), "NA", density_val),
              ifelse(is.na(chao2_val), "NA", chao2_val),
              ifelse(is.na(sac_slope), "NA", sac_slope),
              ifelse(is.na(temporal_dist), "NA", temporal_dist),
              ifelse(is.na(spatial_dist), "NA", spatial_dist)), 
      file = out_file, append = TRUE)
}
# if (file.exists(temp_cube_file)) file.remove(temp_cube_file) # No longer needed
quit(status = 0)
