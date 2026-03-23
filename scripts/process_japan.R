#!/usr/bin/env Rscript
#' Calculate density for Japan only

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
  
  cube <- b3gbi::process_cube(cubepath, separator = ",")
  
  total_occ_map <- b3gbi::compute_indicator_workflow(
    data = cube, type = "total_occ", dim_type = "map",
    shapefile_path = shapefilepath, ne_scale = "large", region = continent,
    include_water = TRUE, shapefile_crs = 4326, ci_type = "none"
  )
  
  area <- total_occ_map$original_bbox %>% sf::st_area() %>% units::set_units("km^2")
  total_occ <- sum(total_occ_map$data$diversity_val, na.rm = TRUE)
  density <- as.numeric(total_occ) / as.numeric(area)
  
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
  
  output_dir <- "output/density_results/asia"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  output_file <- paste0(output_dir, "/", site_id, "_density.RData")
  save(density_obj, file = output_file)
  
  return(density)
}

# Process Japan
continent <- "asia"
country_name <- "Japan"
input_base <- "inst/extdata"
shapefiledir <- file.path(input_base, "ramsar_sites_wkt")
inputdir <- file.path(input_base, paste0("ramsar_site_data_100m_", continent))
output_dir <- "output/density_results/asia"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

sites <- list.files(file.path(inputdir, country_name), pattern = "\\.csv$")
cat("Japan has", length(sites), "sites\n")

processed <- 0
passed <- 0

for (site_file in sites) {
  sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
  site_id <- paste0(country_name, "_", sitename)
  
  if (file.exists(paste0(output_dir, "/", site_id, "_density.RData"))) {
    processed <- processed + 1
    next
  }
  
  result <- tryCatch({
    process_one_site(site_file, country_name, inputdir, shapefiledir, continent)
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(result)) {
    processed <- processed + 1
    if (result >= 0.25) passed <- passed + 1
    cat(".", sep="")
  } else {
    cat("x", sep="")
  }
}

cat("\nDone! Processed:", processed, "Passed:", passed, "\n")
