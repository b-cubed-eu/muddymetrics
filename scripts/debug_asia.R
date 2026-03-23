#!/usr/bin/env Rscript
# Debug version - process Asia with verbose output

library(b3gbi)
library(sf)
library(units)
library(dplyr)
library(purrr)

calculate_density_b3gbi <- function(args_vec) {
  site_file <- args_vec[1]
  country_name <- args_vec[2]
  inputdir <- args_vec[3]
  shapefiledir <- args_vec[4]
  continent <- args_vec[5]
  
  tryCatch({
    sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
    site_id <- paste0(country_name, "_", sitename)
    
    cat("Starting:", site_id, "\n")
    
    shapefilepath <- file.path(shapefiledir, country_name, paste0(sitename, ".wkt"))
    cubepath <- file.path(inputdir, country_name, site_file)
    
    if (!file.exists(shapefilepath) || !file.exists(cubepath)) {
      cat("Missing files for:", site_id, "\n")
      return(NULL)
    }
    
    cube <- b3gbi::process_cube(cubepath, separator = ",")
    cat("Cube processed for:", site_id, "\n")
    
    total_occ_map <- b3gbi::compute_indicator_workflow(
      data = cube, type = "total_occ", dim_type = "map",
      shapefile_path = shapefilepath, ne_scale = "large", region = continent,
      include_water = TRUE, shapefile_crs = 4326, ci_type = "none"
    )
    cat("Indicator computed for:", site_id, "\n")
    
    area <- total_occ_map$original_bbox %>% sf::st_area() %>% units::set_units("km^2")
    total_occ <- sum(total_occ_map$data$diversity_val, na.rm = TRUE)
    density <- as.numeric(total_occ) / as.numeric(area)
    
    density_obj <- list(
      site_id = site_id, sitename = sitename, country = country_name,
      density_km2 = density, total_occurrences = as.numeric(total_occ),
      area_km2 = as.numeric(area), threshold = 0.25, passes_threshold = density >= 0.25,
      timestamp = Sys.time()
    )
    
    output_base <- "output/density_results"
    continent_dir <- file.path(output_base, continent)
    if (!dir.exists(continent_dir)) dir.create(continent_dir, recursive = TRUE)
    
    output_file <- file.path(continent_dir, paste0(site_id, "_density.RData"))
    save(density_obj, file = output_file)
    cat("Saved:", site_id, "\n")
    
  }, error = function(e) {
    cat("Error for", site_file, ":", conditionMessage(e), "\n")
  })
}

# Test with one site
input_base <- "inst/extdata"
shapefiledir <- file.path(input_base, "ramsar_sites_wkt")
inputdir <- file.path(input_base, "ramsar_site_data_100m_asia")
continent_arg <- "asia"

# Find an unprocessed site
existing <- list.files(file.path("output/density_results", continent_arg), pattern = "_density.RData")
processed_ids <- sub("_density.RData$", "", existing)
cat("Already processed:", length(processed_ids), "\n")

countrylist <- list.files(inputdir)
for (country_name in countrylist) {
  country_input_dir <- file.path(inputdir, country_name)
  if (!dir.exists(country_input_dir)) next
  
  sitelist <- list.files(country_input_dir, pattern = ".csv$")
  
  for (site_file in sitelist) {
    sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
    site_id <- paste0(country_name, "_", sitename)
    
    if (site_id %in% processed_ids) next
    
    cat("Found unprocessed:", site_id, "\n")
    args_vec <- c(site_file, country_name, inputdir, shapefiledir, continent_arg)
    result <- calculate_density_b3gbi(args_vec)
    break
  }
  break
}
