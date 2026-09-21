#!/usr/bin/env Rscript
#' Fast density calculation using direct file parsing
#' 
#' This version calculates density directly from CSV files (counting rows for 
#' occurrences) and WKT files (calculating area with sf), avoiding the slow
#' b3gbi processing pipeline.

library(sf)
library(units)
library(dplyr)
library(purrr)

get_density_fast <- function(cubepath, shapefilepath) {
  tryCatch({
    occ_count <- nrow(read.csv(cubepath))
    
    wkt <- readLines(shapefilepath)
    geom <- st_as_sfc(wkt)
    area <- st_area(geom) %>% set_units("km^2")
    
    density <- occ_count / as.numeric(area)
    
    list(
      total_occurrences = occ_count,
      area_km2 = as.numeric(area),
      density_records_km2 = density,
      error = NA_character_
    )
  }, error = function(e) {
    list(
      total_occurrences = NA_real_,
      area_km2 = NA_real_,
      density_records_km2 = NA_real_,
      error = conditionMessage(e)
    )
  })
}

# Actually, let me use a simpler approach: estimate density from the grid cells
# Each CSV has rows representing grid cells with occurrence counts

get_density_from_cube <- function(cubepath, shapefilepath) {
  tryCatch({
    cube <- read.csv(cubepath)
    
    if ("occurrences" %in% names(cube)) {
      total_occ <- sum(cube$occurrences, na.rm = TRUE)
    } else {
      total_occ <- nrow(cube)
    }
    
    wkt <- readLines(shapefilepath)
    geom <- st_as_sfc(wkt)
    area <- st_area(geom) %>% set_units("km^2")
    
    density <- total_occ / as.numeric(area)
    
    list(
      total_occurrences = total_occ,
      area_km2 = as.numeric(area),
      density_records_km2 = density,
      error = NA_character_
    )
  }, error = function(e) {
    list(
      total_occurrences = NA_real_,
      area_km2 = NA_real_,
      density_records_km2 = NA_real_,
      error = conditionMessage(e)
    )
  })
}

process_continent <- function(continent_dir, continent_name, shapefiledir, outputdir) {
  message(paste0("\n===== Processing ", continent_name, " ====="))
  
  countrylist <- list.files(continent_dir)
  all_results <- data.frame()
  
  for (country_name in countrylist) {
    country_input_dir <- file.path(continent_dir, country_name)
    if (!dir.exists(country_input_dir)) next
    
    sitelist <- list.files(country_input_dir, pattern = "\\.csv$")
    
    for (site_file in sitelist) {
      sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
      site_id <- paste0(country_name, "_", sitename)
      
      shapefilepath <- file.path(shapefiledir, country_name, paste0(sitename, ".wkt"))
      cubepath <- file.path(country_input_dir, site_file)
      
      if (!file.exists(shapefilepath)) {
        result <- data.frame(
          site_id = site_id, country = country_name, site_name = sitename,
          total_occurrences = NA, area_km2 = NA, density_records_km2 = NA,
          passes_threshold = NA, error = "missing_wkt", stringsAsFactors = FALSE
        )
      } else {
        d <- get_density_from_cube(cubepath, shapefilepath)
        result <- data.frame(
          site_id = site_id, country = country_name, site_name = sitename,
          total_occurrences = d$total_occurrences, area_km2 = d$area_km2,
          density_records_km2 = d$density_records_km2,
          passes_threshold = d$density_records_km2 >= 0.25,
          error = d$error, stringsAsFactors = FALSE
        )
      }
      
      all_results <- rbind(all_results, result)
    }
    
    message(paste0("Processed ", country_name, " (", length(sitelist), " sites)"))
  }
  
  write.csv(all_results, file = file.path(outputdir, paste0(continent_name, "_density.csv")), row.names = FALSE)
  message(paste0("Saved ", nrow(all_results), " sites to ", continent_name, "_density.csv"))
  
  return(all_results)
}

main <- function() {
  input_base <- "inst/extdata"
  shapefiledir <- file.path(input_base, "ramsar_sites_wkt")
  outputdir <- "output/data_sufficiency"
  
  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE)
  }
  
  continents <- list(
    africa = file.path(input_base, "ramsar_site_data_100m_africa"),
    antarctica = file.path(input_base, "ramsar_site_data_100m_antarctica"),
    asia = file.path(input_base, "ramsar_site_data_100m_asia"),
    europe = file.path(input_base, "ramsar_site_data_100m_europe"),
    northamerica = file.path(input_base, "ramsar_site_data_100m_northamerica"),
    oceania = file.path(input_base, "ramsar_site_data_100m_oceania"),
    southamerica = file.path(input_base, "ramsar_site_data_100m_southamerica")
  )
  
  all_results <- lapply(names(continents), function(cont_name) {
    if (dir.exists(continents[[cont_name]])) {
      process_continent(continents[[cont_name]], cont_name, shapefiledir, outputdir)
    } else {
      NULL
    }
  })
  
  global_results <- do.call(rbind, Filter(Negate(is.null), all_results))
  
  write.csv(global_results, file = file.path(outputdir, "all_sites_density.csv"), row.names = FALSE)
  
  passing_sites <- global_results %>%
    filter(passes_threshold == TRUE) %>%
    arrange(desc(density_records_km2))
  
  write.csv(passing_sites, file = file.path(outputdir, "sites_passing_density_threshold.csv"), row.names = FALSE)
  
  summary_stats <- list(
    total_sites = nrow(global_results),
    passing_sites = nrow(passing_sites),
    failing_sites = nrow(global_results) - nrow(passing_sites),
    pass_rate = nrow(passing_sites) / nrow(global_results),
    threshold = 0.25,
    mean_density = mean(global_results$density_records_km2, na.rm = TRUE),
    median_density = median(global_results$density_records_km2, na.rm = TRUE)
  )
  
  message("\n===== SUMMARY =====")
  message(paste0("Total sites: ", summary_stats$total_sites))
  message(paste0("Passing (≥0.25/km²): ", summary_stats$passing_sites))
  message(paste0("Failing: ", summary_stats$failing_sites))
  message(paste0("Pass rate: ", round(summary_stats$pass_rate * 100, 1), "%"))
  message(paste0("Mean density: ", round(summary_stats$mean_density, 4)))
  message(paste0("Median density: ", round(summary_stats$median_density, 4)))
  
  return(invisible(summary_stats))
}

main()
