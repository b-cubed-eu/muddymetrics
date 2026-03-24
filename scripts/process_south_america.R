library(sf)
library(units)
library(dplyr)

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

continent_dir <- "inst/extdata/ramsar_site_data_100m_southamerica"
continent_name <- "southamerica"
shapefiledir <- "inst/extdata/ramsar_sites_wkt"
outputdir <- "output/data_sufficiency"

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
message(paste0("Saved South America: ", nrow(all_results), " sites"))
