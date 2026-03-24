#' Calculate and Validate Density Indicators for All Ramsar Sites
#'
#' @description
#' This function calculates the occurrence density (records/km²) for all Ramsar sites,
#' saves the results as .RData files, and validates against the Troia 2016 threshold
#' of 0.25 records/km².
#'
#' @param inputdir Character. Path to the base directory containing country-level
#'   subdirectories with GBIF data cubes.
#' @param outputdir Character. Path to save the density .RData files.
#' @param shapefiledir Character. Path to the base directory containing site WKT files.
#' @param threshold Numeric. The density threshold for classification (default 0.25).
#'
#' @return A data frame with site names, density values, and pass/fail status.
#'
#' @importFrom b3gbi process_cube compute_indicator_workflow
#' @importFrom sf st_area st_read
#' @importFrom units set_units
#' @export
calculate_density_indicators <- function(inputdir,
                                        outputdir,
                                        shapefiledir,
                                        threshold = 0.25) {

  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE)
  }

  density_results <- data.frame(
    site_id = character(),
    country = character(),
    density_records_km2 = numeric(),
    passes_threshold = logical(),
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
        cubepath <- file.path(country_input_dir, site_file)

        if (!file.exists(shapefilepath) || !file.exists(cubepath)) {
          message(paste0("Skipping ", site_id, ": missing files"))
          next
        }

        cube <- b3gbi::process_cube(cubepath, separator = ",")

        total_occ_map <- b3gbi::compute_indicator_workflow(
          data = cube,
          type = "total_occ",
          dim_type = "map",
          shapefile_path = shapefilepath,
          ne_scale = "large",
          region = "Africa",
          include_water = TRUE,
          shapefile_crs = 4326,
          ci_type = "none"
        )

        density_val <- get_ramsar_occ_density(total_occ_map)
        density_numeric <- as.numeric(density_val)

        passes <- density_numeric >= threshold

        density_obj <- list(
          site_id = site_id,
          sitename = sitename,
          country = country_name,
          density = density_val,
          density_numeric = density_numeric,
          threshold = threshold,
          passes_threshold = passes,
          timestamp = Sys.time()
        )

        output_file <- file.path(outputdir, paste0(site_id, "_density.RData"))
        save(density_obj, file = output_file)

        density_results <- rbind(density_results, data.frame(
          site_id = site_id,
          country = country_name,
          density_records_km2 = density_numeric,
          passes_threshold = passes,
          stringsAsFactors = FALSE
        ))

      }, error = function(e) {
        message(paste0("Error processing ", site_file, ": ", conditionMessage(e)))
      })
    }
  }

  return(density_results)
}


#' Validate Density Results Against Troia 2016 Threshold
#'
#' @param density_results Data frame with density calculation results.
#' @param threshold Numeric. The density threshold (default 0.25).
#'
#' @return A list with summary statistics and classification.
#'
#' @export
validate_density_results <- function(density_results, threshold = 0.25) {

  total_sites <- nrow(density_results)
  passing_sites <- sum(density_results$passes_threshold, na.rm = TRUE)
  failing_sites <- total_sites - passing_sites

  if (total_sites == 0) {
    return(list(
      total_sites = 0,
      passing_sites = 0,
      failing_sites = 0,
      pass_rate = NaN,
      threshold = threshold,
      mean_density = NA_real_,
      median_density = NA_real_,
      min_density = NA_real_,
      max_density = NA_real_
    ))
  }

  summary_stats <- list(
    total_sites = total_sites,
    passing_sites = passing_sites,
    failing_sites = failing_sites,
    pass_rate = passing_sites / total_sites,
    threshold = threshold,
    mean_density = mean(density_results$density_records_km2, na.rm = TRUE),
    median_density = median(density_results$density_records_km2, na.rm = TRUE),
    min_density = min(density_results$density_records_km2, na.rm = TRUE),
    max_density = max(density_results$density_records_km2, na.rm = TRUE)
  )

  return(summary_stats)
}
