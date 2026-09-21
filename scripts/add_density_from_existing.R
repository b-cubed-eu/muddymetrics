#!/usr/bin/env Rscript
#' Add density for sites that were already successfully processed

library(b3gbi)
library(sf)
library(units)

add_density_from_existing <- function(continent) {
  metric_dir <- paste0("output/ramsar_metric_results_100m/", continent)
  density_dir <- paste0("output/density_results/", continent)
  
  if (!dir.exists(density_dir)) dir.create(density_dir, recursive = TRUE)
  
  # Get all country directories
  countries <- list.files(metric_dir)
  
  # Build list of files to process
  file_list <- list()
  idx <- 1
  for (c in countries) {
    country_dir <- paste0(metric_dir, "/", c)
    if (!dir.exists(country_dir)) next
    files <- list.files(country_dir, pattern = "total_occ_map.RData$")
    for (f in files) {
      file_list[[idx]] <- list(country_dir=country_dir, country=c, filename=f)
      idx <- idx + 1
    }
  }
  
  message(paste0("Found ", length(file_list), " existing indicator files for ", continent))
  
  added <- 0
  skipped <- 0
  
  for (i in seq_along(file_list)) {
    item <- file_list[[i]]
    country_dir <- item$country_dir
    country <- item$country
    filename <- item$filename
    
    # Extract site_id from filename
    site_id <- sub("_total_occ_map.RData$", "", filename)
    
    # Check if already have density
    density_file <- paste0(density_dir, "/", site_id, "_density.RData")
    if (file.exists(density_file)) {
      skipped <- skipped + 1
      next
    }
    
    # Full path to the file
    full_path <- paste0(country_dir, "/", filename)
    
    # Try to load the total_occ_map
    result <- tryCatch({
      e <- new.env()
      suppressWarnings(load(full_path, envir = e))
      
      # Get the object
      obj_name <- ls(e)[1]
      if (is.null(obj_name)) return(NULL)
      
      total_occ_map <- get(obj_name, e)
      
      # Calculate density
      area <- total_occ_map$original_bbox %>% sf::st_area() %>% units::set_units("km^2")
      total_occ <- sum(total_occ_map$data$diversity_val, na.rm = TRUE)
      density <- as.numeric(total_occ) / as.numeric(area)
      
      # Create density object
      density_obj <- list(
        site_id = site_id,
        sitename = site_id,
        country = country,
        density_km2 = density,
        total_occurrences = as.numeric(total_occ),
        area_km2 = as.numeric(area),
        threshold = 0.25,
        passes_threshold = density >= 0.25,
        timestamp = Sys.time()
      )
      
      save(density_obj, file = density_file)
      
      return(TRUE)
    }, error = function(e) {
      message(paste0("Error loading ", filename, ": ", conditionMessage(e)))
      return(NULL)
    })
    
    if (!is.null(result)) {
      added <- added + 1
      if (added %% 50 == 0) {
        message(paste0("Added ", added, " density results..."))
      }
    } else {
      skipped <- skipped + 1
    }
  }
  
  message(paste0("Added: ", added, ", Skipped: ", skipped))
}

# Run for Asia
add_density_from_existing("asia")
