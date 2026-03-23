#!/usr/bin/env Rscript
#' Calculate density using b3gbi with parallel processing

library(b3gbi)
library(sf)
library(units)
library(dplyr)
library(purrr)
library(parallel)

calculate_density_b3gbi <- function(args_vec) {
  site_file <- args_vec[1]
  country_name <- args_vec[2]
  inputdir <- args_vec[3]
  shapefiledir <- args_vec[4]
  continent <- args_vec[5]
  
  tryCatch({
    sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
    site_id <- paste0(country_name, "_", sitename)
    
    shapefilepath <- file.path(shapefiledir, country_name, paste0(sitename, ".wkt"))
    cubepath <- file.path(inputdir, country_name, site_file)
    
    if (!file.exists(shapefilepath) || !file.exists(cubepath)) {
      return(data.frame(
        site_id = site_id, country = country_name, site_name = sitename,
        total_occurrences = NA_real_, area_km2 = NA_real_, density_records_km2 = NA_real_,
        passes_threshold = NA, error = "missing_files", stringsAsFactors = FALSE
      ))
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
    
    data.frame(
      site_id = site_id, country = country_name, site_name = sitename,
      total_occurrences = as.numeric(total_occ), area_km2 = as.numeric(area),
      density_records_km2 = density, passes_threshold = density >= 0.25,
      error = NA_character_, stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      site_id = paste0(country_name, "_", sub("_data$", "", tools::file_path_sans_ext(site_file))),
      country = country_name, site_name = sub("_data$", "", tools::file_path_sans_ext(site_file)),
      total_occurrences = NA_real_, area_km2 = NA_real_, density_records_km2 = NA_real_,
      passes_threshold = NA, error = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
}

args <- commandArgs(trailingOnly = TRUE)
continent_arg <- if (length(args) > 0) args[1] else "europe"

input_base <- "inst/extdata"
shapefiledir <- file.path(input_base, "ramsar_sites_wkt")
inputdir <- file.path(input_base, paste0("ramsar_site_data_100m_", continent_arg))
output_base <- "output/density_results"
continent_dir <- file.path(output_base, continent_arg)

existing <- if (dir.exists(continent_dir)) {
  list.files(continent_dir, pattern = "_density\\.RData$")
} else {
  character()
}
processed_ids <- sub("_density\\.RData$", "", existing)
message(paste0("Already processed: ", length(processed_ids), " sites"))

countrylist <- list.files(inputdir)
tasks <- list()

for (country_name in countrylist) {
  country_input_dir <- file.path(inputdir, country_name)
  if (!dir.exists(country_input_dir)) next
  
  sitelist <- list.files(country_input_dir, pattern = "\\.csv$")
  
  for (site_file in sitelist) {
    sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
    site_id <- paste0(country_name, "_", sitename)
    
    if (!site_id %in% processed_ids) {
      tasks[[length(tasks) + 1]] <- c(site_file, country_name, inputdir, shapefiledir, continent_arg)
    }
  }
}

message(paste0("Remaining tasks: ", length(tasks)))

if (length(tasks) == 0) {
  message("Nothing to do!")
  quit()
}

message("Starting parallel processing with 4 cores...")

# Use mclapply for parallel processing
results_list <- mclapply(tasks, calculate_density_b3gbi, mc.cores = 4)

all_results <- do.call(rbind, results_list)

outputdir <- "output/data_sufficiency"
if (!dir.exists(outputdir)) dir.create(outputdir, recursive = TRUE)

write.csv(all_results, file.path(outputdir, paste0(continent_arg, "_density_b3gbi.csv")), 
          row.names = FALSE, append = TRUE)

message(paste0("\nCompleted: ", nrow(all_results), " sites"))
message(paste0("Passing: ", sum(all_results$passes_threshold == TRUE, na.rm=TRUE)))
