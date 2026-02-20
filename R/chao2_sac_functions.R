#' Calculate Chao2 Species Richness Estimator
#'
#' @description
#' Calculates the Chao2 estimator of species richness for incidence data.
#' Chao2 estimates the true species richness based on the number of rare species
#' (species observed in 1 or 2 samples).
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param shapefile_path Optional path to WKT file for spatial filtering.
#'
#' @return A list with:
#'   - chao2: The Chao2 estimator
#'   - observed: Number of observed species
#'   - completeness: Inventory completeness ratio (observed/chao2)
#'   - f1: Number of species in 1 sample
#'   - f2: Number of species in 2 samples
#'
#' @export
calculate_chao2 <- function(cube, shapefile_path = NULL) {

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  species_col <- cube$specieskey

  valid_species <- species_col[!is.na(species_col)]

  if (length(valid_species) == 0) {
    return(list(
      chao2 = NA_real_,
      observed = 0,
      completeness = NA_real_,
      f1 = 0,
      f2 = 0
    ))
  }

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

  completeness <- if (chao2_est > 0) S_obs / chao2_est else NA_real_

  return(list(
    chao2 = chao2_est,
    observed = S_obs,
    completeness = completeness,
    f1 = f1,
    f2 = f2
  ))
}


#' Calculate Species Accumulation Curve (SAC) Slope
#'
#' @description
#' Calculates the slope of the species accumulation curve as a measure of
#' survey saturation. Lower slopes indicate more complete sampling.
#' Target threshold: <= 0.10 (Troia 2016).
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param shapefile_path Optional path to WKT file for spatial filtering.
#' @param n_iterations Number of randomizations for SAC (default 50).
#'
#' @return A list with:
#'   - slope: The SAC slope
#'   - intercept: The SAC intercept
#'   - r_squared: Goodness of fit
#'
#' @export
calculate_sac_slope <- function(cube,
                                 shapefile_path = NULL,
                                 n_iterations = 50) {

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  unique_occurrences <- unique(cube$occurrenceId)

  if (length(unique_occurrences) < 10) {
    return(list(
      slope = NA_real_,
      intercept = NA_real_,
      r_squared = NA_real_
    ))
  }

  occ_matrix <- table(cube$occurrenceId, cube$specieskey)

  if (nrow(occ_matrix) < 10 || ncol(occ_matrix) < 2) {
    return(list(
      slope = NA_real_,
      intercept = NA_real_,
      r_squared = NA_real_
    ))
  }

  sac_result <- vegan::specaccum(
    occ_matrix,
    method = "random",
    permutations = n_iterations
  )

  sites <- sac_result$sites
  richness <- sac_result$richness

  valid_idx <- !is.na(richness) & sites > 0

  if (sum(valid_idx) < 2) {
    return(list(
      slope = NA_real_,
      intercept = NA_real_,
      r_squared = NA_real_
    ))
  }

  lm_result <- lm(richness[valid_idx] ~ sites[valid_idx])

  slope <- coef(lm_result)[2]
  intercept <- coef(lm_result)[1]
  r_squared <- summary(lm_result)$r.squared

  return(list(
    slope = slope,
    intercept = intercept,
    r_squared = r_squared
  ))
}


#' Batch Calculate Inventory Completeness Metrics
#'
#' @param inputdir Character. Path to directory containing site CSV files.
#' @param outputdir Character. Path to save results.
#' @param shapefiledir Character. Path to WKT files.
#'
#' @return A data frame with site IDs and completeness metrics.
#'
#' @export
calculate_inventory_completeness_batch <- function(inputdir,
                                                    outputdir,
                                                    shapefiledir) {

  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE)
  }

  results <- data.frame(
    site_id = character(),
    country = character(),
    chao2 = numeric(),
    observed_species = numeric(),
    completeness_chao2 = numeric(),
    sac_slope = numeric(),
    sac_r2 = numeric(),
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
          next
        }

        cubepath <- file.path(country_input_dir, site_file)

        chao2_result <- calculate_chao2(cubepath, shapefilepath)
        sac_result <- calculate_sac_slope(cubepath, shapefilepath)

        results <- rbind(results, data.frame(
          site_id = site_id,
          country = country_name,
          chao2 = chao2_result$chao2,
          observed_species = chao2_result$observed,
          completeness_chao2 = chao2_result$completeness,
          sac_slope = sac_result$slope,
          sac_r2 = sac_result$r_squared,
          stringsAsFactors = FALSE
        ))

      }, error = function(e) {
        message(paste0("Error processing ", site_file, ": ", conditionMessage(e)))
      })
    }
  }

  return(results)
}
