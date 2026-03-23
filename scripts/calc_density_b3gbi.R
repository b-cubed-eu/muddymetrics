#!/usr/bin/env Rscript
#' Calculate density - process one site at a time, save immediately

library(b3gbi)
library(sf)
library(units)

process_one_site <- function(site_file, country_name, inputdir, shapefiledir, continent) {
  sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
  site_id <- paste0(country_name, "_", sitename)
  
  shapefilepath <- file.path(shapefiledir, country_name, paste0(sitename, ".wkt"))
  cubepath <- file.path(inputdir, country_name, site_file)
  
  if (!file.exists(shapefilepath) || !file.exists(cubepath)) {
    return(NULL)
  }
  
  # Process one site
  cube <- b3gbi::process_cube(cubepath, separator = ",")
  
  total_occ_map <- b3gbi::compute_indicator_workflow(
    data = cube, type = "total_occ", dim_type = "map",
    shapefile_path = shapefilepath, ne_scale = "large", region = continent,
    include_water = TRUE, shapefile_crs = 4326, ci_type = "none"
  )
  
  area <- total_occ_map$original_bbox %>% sf::st_area() %>% units::set_units("km^2")
  total_occ <- sum(total_occ_map$data$diversity_val, na.rm = TRUE)
  density <- as.numeric(total_occ) / as.numeric(area)
  
  # Save immediately
  density_obj <- list(
    site_id = site_id,
    sitename = sitename,
    country = country_name,
    density_km2 = density,
    total_occurrences = as.numeric(total_occ),
    area_km2 = as.numeric(area),
    threshold = 0.25,
    passes_threshold = density >= 0.25,
    timestamp = Sys.time()
  )
  
  output_dir <- paste0("output/density_results/", continent)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_file <- paste0(output_dir, "/", site_id, "_density.RData")
  save(density_obj, file = output_file)
  
  return(density)
}

# Main processing
args <- commandArgs(trailingOnly = TRUE)
continent <- if (length(args) > 0) args[1] else "asia"

input_base <- "inst/extdata"
shapefiledir <- file.path(input_base, "ramsar_sites_wkt")
inputdir <- file.path(input_base, paste0("ramsar_site_data_100m_", continent))
output_dir <- paste0("output/density_results/", continent)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Get all countries
countries <- sort(list.files(inputdir))
total_processed <- 0
total_passed <- 0
errors <- 0

for (country_name in countries) {
  country_dir <- file.path(inputdir, country_name)
  if (!dir.exists(country_dir)) next
  
  # Check how many already processed in this country
  existing <- list.files(output_dir, pattern = paste0("^", country_name, "_.*_density.RData$"))
  
  sites <- list.files(country_dir, pattern = "\\.csv$")
  to_process <- sites
  
  if (length(existing) >= length(sites)) {
    cat(country_name, ": already done (", length(existing), "/", length(sites), ")\n", sep="")
    next
  }
  
  cat("Processing ", country_name, " (", length(sites), " sites)...\n", sep="")
  
  for (site_file in sites) {
    sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
    site_id <- paste0(country_name, "_", sitename)
    
    # Check if already processed
    if (file.exists(paste0(output_dir, "/", site_id, "_density.RData"))) {
      total_processed <- total_processed + 1
      next
    }
    
    # Try to process
    result <- tryCatch({
      process_one_site(site_file, country_name, inputdir, shapefiledir, continent)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(result)) {
      total_processed <- total_processed + 1
      if (result >= 0.25) total_passed <- total_passed + 1
    } else {
      errors <- errors + 1
    }
  }
  
  cat("  Done ", country_name, " (processed:", total_processed, ", passed:", total_passed, ")\n", sep="")
}

cat("\n===== COMPLETE =====\n")
cat("Total processed:", total_processed, "\n")
cat("Total passed:", total_passed, "\n")
cat("Errors:", errors, "\n")
