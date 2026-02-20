#' Calculate Alpha Diversity Metrics
#'
#' @description
#' Calculates various alpha diversity metrics including Hill numbers (0, 1, 2),
#' Shannon diversity, and Simpson diversity for a Ramsar site.
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param shapefile_path Optional path to WKT file for spatial filtering.
#'
#' @return A list with:
#'   - hill_0: Species richness (Hill number 0)
#'   - hill_1: Shannon diversity (Hill number 1)
#'   - hill_2: Simpson diversity (Hill number 2)
#'   - shannon_h: Shannon entropy
#'   - simpson_d: Simpson diversity index
#'   - simpson_1d: Simpson reciprocal index
#'
#' @export
calculate_alpha_diversity <- function(cube, shapefile_path = NULL) {

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  species_col <- cube$specieskey
  occurrences_col <- cube$occurrences

  valid_data <- !is.na(species_col) & !is.na(occurrences_col) & occurrences_col > 0

  if (!any(valid_data)) {
    return(list(
      hill_0 = NA_integer_,
      hill_1 = NA_real_,
      hill_2 = NA_real_,
      shannon_h = NA_real_,
      simpson_d = NA_real_,
      simpson_1d = NA_real_
    ))
  }

  species_counts <- tapply(
    occurrences_col[valid_data],
    species_col[valid_data],
    sum
  )

  species_counts <- species_counts[species_counts > 0]

  total_abundance <- sum(species_counts)

  prop_abundances <- species_counts / total_abundance

  hill_0 <- length(species_counts)

  if (hill_0 == 0) {
    return(list(
      hill_0 = 0,
      hill_1 = 0,
      hill_2 = 0,
      shannon_h = 0,
      simpson_d = 0,
      simpson_1d = 0
    ))
  }

  shannon_h <- -sum(prop_abundances * log(prop_abundances))

  hill_1 <- exp(shannon_h)

  simpson_d <- 1 - sum(prop_abundances^2)

  simpson_1d <- 1 / sum(prop_abundances^2)

  lambda <- sum(prop_abundances^2)
  hill_2 <- 1 / lambda

  return(list(
    hill_0 = hill_0,
    hill_1 = hill_1,
    hill_2 = hill_2,
    shannon_h = shannon_h,
    simpson_d = simpson_d,
    simpson_1d = simpson_1d
  ))
}


#' Calculate Evenness Metrics
#'
#' @description
#' Calculates Pielou's Evenness and Williams' Evenness for a Ramsar site.
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param shapefile_path Optional path to WKT file for spatial filtering.
#'
#' @return A list with:
#'   - pielou_j: Pielou's evenness (J)
#'   - williams_w: Williams' evenness (W)
#'
#' @export
calculate_evenness <- function(cube, shapefile_path = NULL) {

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  species_col <- cube$specieskey
  occurrences_col <- cube$occurrences

  valid_data <- !is.na(species_col) & !is.na(occurrences_col) & occurrences_col > 0

  if (!any(valid_data)) {
    return(list(
      pielou_j = NA_real_,
      williams_w = NA_real_
    ))
  }

  species_counts <- tapply(
    occurrences_col[valid_data],
    species_col[valid_data],
    sum
  )

  species_counts <- species_counts[species_counts > 0]

  total_abundance <- sum(species_counts)
  s <- length(species_counts)

  if (s <= 1 || total_abundance <= 1) {
    return(list(
      pielou_j = NA_real_,
      williams_w = NA_real_
    ))
  }

  prop_abundances <- species_counts / total_abundance

  shannon_h <- -sum(prop_abundances * log(prop_abundances))

  pielou_j <- shannon_h / log(s)

  var_s <- var(species_counts)
  williams_w <- (shannon_h - log(s)) / (s * log(total_abundance))

  return(list(
    pielou_j = pielou_j,
    williams_w = williams_w
  ))
}


#' Calculate Rarity Metrics
#'
#' @description
#' Calculates Area-Based Rarity and Abundance-Based Rarity metrics.
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param site_area_km2 Numeric. Area of the Ramsar site in km².
#' @param shapefile_path Optional path to WKT file for spatial filtering.
#'
#' @return A list with:
#'   - area_rarity: Rarity score based on restricted range
#'   - abundance_rarity: Rarity based on low abundance
#'
#' @export
calculate_rarity_metrics <- function(cube,
                                     site_area_km2 = NULL,
                                     shapefile_path = NULL) {

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  species_col <- cube$specieskey
  occurrences_col <- cube$occurrences
  geometry_col <- cube$geometry

  valid_data <- !is.na(species_col) & !is.na(occurrences_col) & occurrences_col > 0

  if (!any(valid_data)) {
    return(list(
      area_rarity = NA_real_,
      abundance_rarity = NA_real_,
      mean_range_size_km2 = NA_real_
    ))
  }

  species_data <- data.frame(
    specieskey = species_col[valid_data],
    occurrences = occurrences_col[valid_data],
    geometry = geometry_col[valid_data]
  )

  species_abundances <- tapply(
    species_data$occurrences,
    species_data$specieskey,
    sum
  )

  total_occurrences <- sum(species_abundances)

  abundance_rarity <- sum(species_abundances < (total_occurrences * 0.01)) / length(species_abundances)

  unique_cells <- length(unique(species_data$geometry))

  if (!is.null(site_area_km2) && site_area_km2 > 0) {
    area_per_cell <- site_area_km2 / unique_cells
  } else {
    area_per_cell <- NA_real_
  }

  mean_range <- unique_cells * area_per_cell

  species_cells <- tapply(
    rep(1, nrow(species_data)),
    species_data$specieskey,
    sum
  )

  rarity_contribution <- 1 / (species_cells + 1)
  area_rarity <- sum(rarity_contribution) / length(rarity_contribution)

  return(list(
    area_rarity = area_rarity,
    abundance_rarity = abundance_rarity,
    mean_range_size_km2 = mean_range
  ))
}


#' Batch Calculate All Ecological Metrics
#'
#' @param inputdir Character. Path to directory containing site CSV files.
#' @param outputdir Character. Path to save results.
#' @param shapefiledir Character. Path to WKT files.
#'
#' @return A data frame with all ecological metrics per site.
#'
#' @export
calculate_ecological_metrics_batch <- function(inputdir,
                                               outputdir,
                                               shapefiledir) {

  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE)
  }

  results <- data.frame(
    site_id = character(),
    country = character(),
    hill_0 = numeric(),
    hill_1 = numeric(),
    hill_2 = numeric(),
    shannon_h = numeric(),
    simpson_d = numeric(),
    pielou_j = numeric(),
    abundance_rarity = numeric(),
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

        alpha_div <- calculate_alpha_diversity(cubepath, shapefilepath)
        evenness <- calculate_evenness(cubepath, shapefilepath)
        rarity <- calculate_rarity_metrics(cubepath, shapefilepath = shapefilepath)

        results <- rbind(results, data.frame(
          site_id = site_id,
          country = country_name,
          hill_0 = alpha_div$hill_0,
          hill_1 = alpha_div$hill_1,
          hill_2 = alpha_div$hill_2,
          shannon_h = alpha_div$shannon_h,
          simpson_d = alpha_div$simpson_d,
          pielou_j = evenness$pielou_j,
          abundance_rarity = rarity$abundance_rarity,
          stringsAsFactors = FALSE
        ))

      }, error = function(e) {
        message(paste0("Error processing ", site_file, ": ", conditionMessage(e)))
      })
    }
  }

  return(results)
}
