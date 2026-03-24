#' Calculate Phylogenetic Diversity
#'
#' @description
#' Wrapper function to calculate phylogenetic diversity indicators using pdindicatoR.
#' Requires a phylogenetic tree (Newick or Nexus format) and species occurrence data.
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param tree_path Character. Path to phylogenetic tree file (.nex or .newick).
#' @param shapefile_path Optional path to WKT file for spatial filtering.
#'
#' @return A list with phylogenetic diversity metrics.
#'
#' @export
calculate_phylogenetic_diversity <- function(cube,
                                            tree_path,
                                            shapefile_path = NULL) {

  if (!requireNamespace("pdindicatoR", quietly = TRUE)) {
    stop("pdindicatoR package required. Install with: install.packages('pdindicatoR', repos = 'https://b-cubed-eu.r-universe.dev')")
  }

  if (!file.exists(tree_path)) {
    stop("Phylogenetic tree file not found: ", tree_path)
  }

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  tree <- ape::read.tree(tree_path)

  species_list <- unique(cube$specieskey)

  species_in_tree <- species_list[species_list %in% tree$tip.label]

  if (length(species_in_tree) < 2) {
    return(list(
      faith_pd = NA_real_,
      pd_coverage = NA_real_,
      species_matched = 0,
      total_species = length(species_list)
    ))
  }

  pd_values <- pdindicatoR::calculate_faithpd(
    tree = tree,
    species = species_in_tree,
    mrca_node_id = NULL
  )

  pd_coverage <- length(species_in_tree) / length(species_list)

  return(list(
    faith_pd = pd_values,
    pd_coverage = pd_coverage,
    species_matched = length(species_in_tree),
    total_species = length(species_list)
  ))
}


#' Calculate Invasive Species Impact
#'
#' @description
#' Wrapper function to calculate invasive species impact indicators using impIndicator.
#' Requires EICAT (Environmental Impact Classification for Alien Taxa) assessment data.
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param eicat_data A data frame with EICAT assessments (species, impact_category, mechanism).
#' @param shapefile_path Optional path to WKT file for spatial filtering.
#' @param method Character. Method for aggregation (default "mean_cumulative").
#'
#' @return A list with invasive impact metrics.
#'
#' @export
calculate_invasive_impact <- function(cube,
                                     eicat_data,
                                     shapefile_path = NULL,
                                     method = "mean_cumulative") {

  if (!requireNamespace("impIndicator", quietly = TRUE)) {
    stop("impIndicator package required. Install with: install.packages('impIndicator', repos = 'https://b-cubed-eu.r-universe.dev')")
  }

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  if (!"species" %in% names(eicat_data)) {
    stop("EICAT data must have a 'species' column")
  }

  species_col <- cube$species

  if (is.null(species_col)) {
    species_col <- cube$specieskey
  }

  cube$species_name <- species_col

  unique_species <- unique(cube$species_name)

  eicat_matched <- eicat_data[eicat_data$species %in% unique_species, ]

  if (nrow(eicat_matched) == 0) {
    return(list(
      total_invasive_species = 0,
      impact_indicator = NA_real_,
      species_assessed = 0,
      total_species = length(unique_species)
    ))
  }

  impact_result <- tryCatch({
    impIndicator::compute_impact_indicator(
      cube = cube,
      impact_data = eicat_matched,
      method = method
    )
  }, error = function(e) {
    message("Error computing impact indicator: ", conditionMessage(e))
    NULL
  })

  if (is.null(impact_result)) {
    return(list(
      total_invasive_species = nrow(eicat_matched),
      impact_indicator = NA_real_,
      species_assessed = nrow(eicat_matched),
      total_species = length(unique_species)
    ))
  }

  return(list(
    total_invasive_species = nrow(eicat_matched),
    impact_indicator = impact_result,
    species_assessed = nrow(eicat_matched),
    total_species = length(unique_species)
  ))
}


#' Batch Calculate Specialized Indicators
#'
#' @param inputdir Character. Path to directory containing site CSV files.
#' @param outputdir Character. Path to save results.
#' @param shapefiledir Character. Path to WKT files.
#' @param tree_path Character. Path to phylogenetic tree file.
#' @param eicat_data Data frame. EICAT assessment data.
#'
#' @return A data frame with specialized indicators.
#'
#' @export
calculate_specialized_batch <- function(inputdir,
                                        outputdir,
                                        shapefiledir,
                                        tree_path = NULL,
                                        eicat_data = NULL) {

  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE)
  }

  results <- data.frame(
    site_id = character(),
    country = character(),
    faith_pd = numeric(),
    pd_coverage = numeric(),
    invasive_species = numeric(),
    impact_indicator = numeric(),
    stringsAsFactors = FALSE
  )

  if (!is.null(tree_path) && !file.exists(tree_path)) {
    message("Warning: Tree path not found. Phylogenetic diversity will be skipped.")
    tree_path <- NULL
  }

  if (!is.null(eicat_data) && !is.data.frame(eicat_data)) {
    message("Warning: EICAT data invalid. Invasive impact will be skipped.")
    eicat_data <- NULL
  }

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

        faith_pd <- NA_real_
        pd_coverage <- NA_real_
        invasive_species <- NA_real_
        impact_indicator <- NA_real_

        if (!is.null(tree_path)) {
          pd_result <- tryCatch({
            calculate_phylogenetic_diversity(cubepath, tree_path, shapefilepath)
          }, error = function(e) {
            message("PD error for ", site_id, ": ", conditionMessage(e))
            NULL
          })

          if (!is.null(pd_result)) {
            faith_pd <- pd_result$faith_pd
            pd_coverage <- pd_result$pd_coverage
          }
        }

        if (!is.null(eicat_data)) {
          inv_result <- tryCatch({
            calculate_invasive_impact(cubepath, eicat_data, shapefilepath)
          }, error = function(e) {
            message("Invasive error for ", site_id, ": ", conditionMessage(e))
            NULL
          })

          if (!is.null(inv_result)) {
            invasive_species <- inv_result$total_invasive_species
            impact_indicator <- inv_result$impact_indicator
          }
        }

        results <- rbind(results, data.frame(
          site_id = site_id,
          country = country_name,
          faith_pd = faith_pd,
          pd_coverage = pd_coverage,
          invasive_species = invasive_species,
          impact_indicator = impact_indicator,
          stringsAsFactors = FALSE
        ))

      }, error = function(e) {
        message(paste0("Error processing ", site_file, ": ", conditionMessage(e)))
      })
    }
  }

  return(results)
}
