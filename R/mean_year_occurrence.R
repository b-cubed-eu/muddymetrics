#' Calculate Mean Year of Occurrence for a Ramsar Site
#'
#' @description
#' Calculates the weighted mean year of occurrence as a proxy for data recency.
#' Higher mean years indicate more recent data.
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param shapefile_path Optional path to the site's WKT/shapefile for filtering.
#'
#' @return A numeric value representing the mean year of occurrence.
#'
#' @export
calculate_mean_year <- function(cube, shapefile_path = NULL) {

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  year_col <- cube$year
  occurrences_col <- cube$occurrences

  valid_data <- !is.na(year_col) & !is.na(occurrences_col) & year_col > 0

  if (!any(valid_data)) {
    return(NA_real_)
  }

  mean_year <- sum(year_col[valid_data] * occurrences_col[valid_data], na.rm = TRUE) /
    sum(occurrences_col[valid_data], na.rm = TRUE)

  return(mean_year)
}


#' Calculate Mean Year for All Sites in a Directory
#'
#' @param inputdir Character. Path to directory containing site CSV files.
#' @param outputdir Character. Path to save results.
#' @param shapefiledir Character. Path to WKT files.
#'
#' @return A data frame with site IDs and mean years.
#'
#' @export
calculate_mean_year_batch <- function(inputdir, outputdir, shapefiledir) {

  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE)
  }

  results <- data.frame(
    site_id = character(),
    country = character(),
    mean_year = numeric(),
    year_min = numeric(),
    year_max = numeric(),
    total_occurrences = numeric(),
    stringsAsFactors = FALSE
  )

  countrylist <- list.files(inputdir)

  for (country_name in countrylist) {
    country_input_dir <- file.path(inputdir, country_name)

    if (!dir.exists(country_input_dir)) next

    sitelist <- list.files(country_input_dir, pattern = "\\.csv$")

    for (site_file in sitelist) {
      tryCatch({
        sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
        site_id <- paste0(country_name, "_", sitename)

        shapefilepath <- file.path(shapefiledir, country_name, paste0(sitename, ".wkt"))

        if (!file.exists(shapefilepath)) {
          message(paste0("Skipping ", site_id, ": missing shapefile"))
          next
        }

        mean_yr <- calculate_mean_year(
          cube = file.path(country_input_dir, site_file),
          shapefile_path = shapefilepath
        )

        cube <- readr::read_csv(
          file.path(country_input_dir, site_file),
          show_col_types = FALSE
        )

        year_min <- min(cube$year, na.rm = TRUE)
        year_max <- max(cube$year, na.rm = TRUE)
        total_occ <- sum(cube$occurrences, na.rm = TRUE)

        results <- rbind(results, data.frame(
          site_id = site_id,
          country = country_name,
          mean_year = mean_yr,
          year_min = year_min,
          year_max = year_max,
          total_occurrences = total_occ,
          stringsAsFactors = FALSE
        ))

      }, error = function(e) {
        message(paste0("Error processing ", site_file, ": ", conditionMessage(e)))
      })
    }
  }

  return(results)
}


#' Filter cube data by shapefile boundary
#'
#' @param cube Data frame. The GBIF occurrence cube.
#' @param shapefile_path Character. Path to WKT file.
#'
#' @return Filtered data frame.
#'
filter_by_shapefile <- function(cube, shapefile_path) {

  site_wkt <- readr::read_file(shapefile_path)

  site_geom <- sf::st_as_sfc(site_wkt)

  cube_coords <- cube$geometry

  coords_matrix <- do.call(rbind, lapply(strsplit(cube_coords, "\\|"), as.numeric))

  colnames(coords_matrix) <- c("lon", "lat")

  pts <- sf::st_as_sf(
    as.data.frame(coords_matrix),
    coords = c("lon", "lat"),
    crs = 4326
  )

  pts_in <- sf::st_intersects(pts, site_geom, sparse = FALSE)

  cube[pts_in, ]
}
