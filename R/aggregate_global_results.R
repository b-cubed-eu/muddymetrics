#' Aggregate Global Results from All Ramsar Sites
#'
#' @description
#' Consolidates all site-level indicator files (.RData) into a single master
#' dataset with all calculated metrics for global analysis.
#'
#' @param input_base Character. Base path to indicator output directories.
#' @param output_dir Character. Path to save the aggregated results.
#'
#' @return A data frame with all site metrics.
#'
#' @export
aggregate_global_results <- function(input_base, output_dir) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  continents <- c("africa", "antarctica", "asia", "europe", "northamerica", "oceania", "southamerica")

  all_sites <- data.frame(
    site_id = character(),
    site_name = character(),
    country = character(),
    continent = character(),
    obs_richness = numeric(),
    total_occurrences = numeric(),
    cum_richness = numeric(),
    density_km2 = numeric(),
    mean_year = numeric(),
    year_min = numeric(),
    year_max = numeric(),
    chao2 = numeric(),
    chao2_completeness = numeric(),
    sac_slope = numeric(),
    stringsAsFactors = FALSE
  )

  for (continent_name in continents) {
    continent_dir <- file.path(input_base, continent_name)

    if (!dir.exists(continent_dir)) {
      message(paste0("Skipping missing continent: ", continent_name))
      next
    }

    country_dirs <- list.files(continent_dir)

    for (country_name in country_dirs) {
      country_dir <- file.path(continent_dir, country_name)

      if (!dir.exists(country_dir)) next

      site_files <- list.files(country_dir, pattern = "\\.RData$")

      for (site_file in site_files) {
        tryCatch({
          site_path <- file.path(country_dir, site_file)

          load_result <- try(load(site_path), silent = TRUE)

          if (inherits(load_result, "try-error")) {
            message(paste0("Could not load: ", site_file))
            next
          }

          site_name <- sub("\\.RData$", "", sub(".*_", "", site_file))

          site_id <- sub("_.*", "", sub("site_", "", site_file))

          new_row <- data.frame(
            site_id = paste0("site_", site_id),
            site_name = site_name,
            country = country_name,
            continent = continent_name,
            obs_richness = NA,
            total_occurrences = NA,
            cum_richness = NA,
            density_km2 = NA,
            mean_year = NA,
            year_min = NA,
            year_max = NA,
            chao2 = NA,
            chao2_completeness = NA,
            sac_slope = NA,
            stringsAsFactors = FALSE
          )

          obj_name <- ls()[sapply(ls(), function(x) {
            obj <- get(x, inherits = FALSE)
            is.list(obj) || is.data.frame(obj)
          })]

          for (obj in obj_name) {
            obj_data <- get(obj)

            if (is.data.frame(obj_data) && "diversity_val" %in% names(obj_data)) {
              new_row$obs_richness <- mean(obj_data$diversity_val, na.rm = TRUE)
            }

            if (is.list(obj_data) && "density_numeric" %in% names(obj_data)) {
              new_row$density_km2 <- obj_data$density_numeric
            }
          }

          all_sites <- rbind(all_sites, new_row)

        }, error = function(e) {
          message(paste0("Error processing ", site_file, ": ", conditionMessage(e)))
        })
      }
    }
  }

  return(all_sites)
}


#' Generate Global Sufficiency Summary CSV
#'
#' @param aggregated_data Data frame from aggregate_global_results().
#' @param output_path Character. Path to save the CSV file.
#' @param threshold_density Numeric. Density threshold (default 0.25).
#' @param threshold_chao2 Numeric. Chao2 completeness threshold (default 0.70).
#' @param threshold_slope Numeric. SAC slope threshold (default 0.10).
#'
#' @return The classified data frame.
#'
#' @export
generate_global_summary <- function(aggregated_data,
                                    output_path,
                                    threshold_density = 0.25,
                                    threshold_chao2 = 0.70,
                                    threshold_slope = 0.10) {

  result <- aggregated_data

  result$data_class <- "Data-Poor"

  sufficient_density <- !is.na(result$density_km2) & result$density_km2 >= threshold_density
  sufficient_chao2 <- !is.na(result$chao2_completeness) & result$chao2_completeness >= threshold_chao2
  sufficient_slope <- !is.na(result$sac_slope) & result$sac_slope <= threshold_slope

  result$data_class[sufficient_density] <- "Data-Rich"

  result$passes_density <- sufficient_density
  result$passes_chao2 <- sufficient_chao2
  result$passes_slope <- sufficient_slope

  result$passes_all_thresholds <- sufficient_density & sufficient_chao2 & sufficient_slope

  result$troia_moderate <- sufficient_density

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(result, output_path)

  message(paste0("Global summary saved to: ", output_path))
  message(paste0("Total sites: ", nrow(result)))
  message(paste0("Data-Rich sites: ", sum(result$data_class == "Data-Rich")))
  message(paste0("Data-Poor sites: ", sum(result$data_class == "Data-Poor")))
  message(paste0("Troia Moderate Threshold (density >= ", threshold_density, "): ",
                 sum(result$troia_moderate)))

  return(result)
}


#' Perform Data Gap Analysis
#'
#' @param summary_data Data frame from generate_global_summary().
#'
#' @return A list with analysis results.
#'
#' @export
perform_data_gap_analysis <- function(summary_data) {

  total_sites <- nrow(summary_data)

  data_rich <- sum(summary_data$data_class == "Data-Rich", na.rm = TRUE)
  data_poor <- sum(summary_data$data_class == "Data-Poor", na.rm = TRUE)

  by_continent <- summary_data |>
    dplyr::group_by(.data$continent) |>
    dplyr::summarise(
      total = dplyr::n(),
      data_rich = sum(.data$data_class == "Data-Rich", na.rm = TRUE),
      data_poor = sum(.data$data_class == "Data-Poor", na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(pct_data_rich = data_rich / .data$total * 100)

  by_country <- summary_data |>
    dplyr::group_by(.data$country) |>
    dplyr::summarise(
      total = dplyr::n(),
      data_rich = sum(.data$data_class == "Data-Rich", na.rm = TRUE),
      data_poor = sum(.data$data_class == "Data-Poor", na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(pct_data_rich = data_rich / .data$total * 100) |>
    dplyr::arrange(dplyr::desc(.data$pct_data_rich))

  threshold_pass <- summary_data |>
    dplyr::summarise(
      density_pass = sum(.data$passes_density, na.rm = TRUE),
      chao2_pass = sum(.data$passes_chao2, na.rm = TRUE),
      slope_pass = sum(.data$passes_slope, na.rm = TRUE),
      all_pass = sum(.data$passes_all_thresholds, na.rm = TRUE)
    )

  return(list(
    total_sites = total_sites,
    data_rich_sites = data_rich,
    data_poor_sites = data_poor,
    pass_rate = data_rich / total_sites,
    by_continent = by_continent,
    by_country = by_country,
    threshold_pass = threshold_pass
  ))
}
