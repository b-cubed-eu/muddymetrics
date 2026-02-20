#' Run European Deep Dive Analysis
#'
#' @description
#' This script runs the complete advanced ecological analysis workflow
#' for European Ramsar sites. It filters for European sites, calculates
#' all advanced metrics, and generates comparative visualizations.
#'
#' @param data_dir Character. Path to GBIF data cubes.
#' @param shapefile_dir Character. Path to WKT boundary files.
#' @param output_dir Character. Path for output files.
#' @param tree_path Character. Optional path to phylogenetic tree.
#' @param eicat_data Data frame. Optional EICAT assessment data.
#'
#' @export
run_european_deep_dive <- function(data_dir,
                                   shapefile_dir,
                                   output_dir,
                                   tree_path = NULL,
                                   eicat_data = NULL) {

  message("Starting European Deep Dive Analysis...")
  message(paste0("Data directory: ", data_dir))
  message(paste0("Output directory: ", output_dir))

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  europe_data_dir <- file.path(data_dir, "europe")

  if (!dir.exists(europe_data_dir)) {
    stop("European data directory not found: ", europe_data_dir)
  }

  message("\n--- Step 1: Calculate Alpha Diversity ---")
  alpha_div_results <- calculate_ecological_metrics_batch(
    inputdir = europe_data_dir,
    outputdir = file.path(output_dir, "alpha_diversity"),
    shapefiledir = shapefile_dir
  )

  readr::write_csv(
    alpha_div_results,
    file.path(output_dir, "europe_alpha_diversity.csv")
  )
  message(paste0("Alpha diversity saved: ", nrow(alpha_div_results), " sites"))

  message("\n--- Step 2: Calculate Inventory Completeness ---")
  inventory_results <- calculate_inventory_completeness_batch(
    inputdir = europe_data_dir,
    outputdir = file.path(output_dir, "inventory"),
    shapefiledir = shapefile_dir
  )

  readr::write_csv(
    inventory_results,
    file.path(output_dir, "europe_inventory_completeness.csv")
  )
  message(paste0("Inventory completeness saved: ", nrow(inventory_results), " sites"))

  message("\n--- Step 3: Calculate Mean Year ---")
  mean_year_results <- calculate_mean_year_batch(
    inputdir = europe_data_dir,
    outputdir = file.path(output_dir, "mean_year"),
    shapefiledir = shapefile_dir
  )

  readr::write_csv(
    mean_year_results,
    file.path(output_dir, "europe_mean_year.csv")
  )
  message(paste0("Mean year saved: ", nrow(mean_year_results), " sites"))

  message("\n--- Step 4: Merge All Results ---")
  merged_results <- alpha_div_results |>
    dplyr::left_join(
      inventory_results |>
        dplyr::select(site_id, chao2, completeness_chao2, sac_slope),
      by = "site_id"
    ) |>
    dplyr::left_join(
      mean_year_results |>
        dplyr::select(site_id, mean_year),
      by = "site_id"
    )

  merged_results <- merged_results |>
    dplyr::mutate(
      data_sufficient = .data$completeness_chao2 >= 0.7 & .data$sac_slope <= 0.10
    )

  readr::write_csv(
    merged_results,
    file.path(output_dir, "europe_complete_analysis.csv")
  )
  message(paste0("Complete analysis saved: ", nrow(merged_results), " sites"))

  message("\n--- Step 5: Generate Comparative Plots ---")

  country_summary <- merged_results |>
    dplyr::group_by(country) |>
    dplyr::summarise(
      n_sites = dplyr::n(),
      mean_hill_1 = mean(.data$hill_1, na.rm = TRUE),
      mean_shannon = mean(.data$shannon_h, na.rm = TRUE),
      mean_pielou = mean(.data$pielou_j, na.rm = TRUE),
      mean_completeness = mean(.data$completeness_chao2, na.rm = TRUE),
      pct_data_sufficient = mean(.data$data_sufficient, na.rm = TRUE) * 100,
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(.data$mean_shannon))

  readr::write_csv(
    country_summary,
    file.path(output_dir, "europe_country_summary.csv")
  )

  if (nrow(merged_results) > 0) {
    plot_diversity_by_country(
      merged_results,
      file.path(output_dir, "europe_diversity_violin.png")
    )

    plot_completeness_by_country(
      merged_results,
      file.path(output_dir, "europe_completeness_violin.png")
    )
  }

  message("\n--- Summary ---")
  message(paste0("Total European sites analyzed: ", nrow(merged_results)))
  message(paste0("Countries represented: ", length(unique(merged_results$country))))
  message(paste0("Data-sufficient sites: ", sum(merged_results$data_sufficient, na.rm = TRUE)))

  message("\nEuropean Deep Dive Analysis Complete!")
  message(paste0("Results saved to: ", output_dir))

  return(merged_results)
}


#' Plot Diversity Metrics by Country
#'
#' @param data Data frame with diversity metrics.
#' @param output_path Character. Path to save plot.
#'
plot_diversity_by_country <- function(data, output_path) {

  data_plot <- data |>
    dplyr::filter(!is.na(.data$hill_1)) |>
    dplyr::mutate(country = forcats::fct_reorder(.data$country, .data$hill_1, .median, .desc = TRUE))

  p <- ggplot2::ggplot(data_plot, ggplot2::aes(x = .data$country, y = .data$hill_1)) |>
    ggplot2::geom_violin(fill = "steelblue", alpha = 0.7) |>
    ggplot2::geom_boxplot(width = 0.2, fill = "white", alpha = 0.5) |>
    ggplot2::theme_minimal() |>
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) |>
    ggplot2::labs(
      title = "Hill Number 1 (Exp Shannon) by Country",
      x = "Country",
      y = "Hill 1 Diversity"
    )

  ggplot2::ggsave(output_path, p, width = 12, height = 6, dpi = 150)
}


#' Plot Inventory Completeness by Country
#'
#' @param data Data frame with completeness metrics.
#' @param output_path Character. Path to save plot.
#'
plot_completeness_by_country <- function(data, output_path) {

  data_plot <- data |>
    dplyr::filter(!is.na(.data$completeness_chao2)) |>
    dplyr::mutate(
      country = forcats::fct_reorder(.data$country, .data$completeness_chao2, .median, .desc = TRUE)
    )

  p <- ggplot2::ggplot(data_plot, ggplot2::aes(x = .data$country, y = .data$completeness_chao2)) |>
    ggplot2::geom_violin(fill = "forestgreen", alpha = 0.7) |>
    ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed", color = "red", linewidth = 1) |>
    ggplot2::geom_boxplot(width = 0.2, fill = "white", alpha = 0.5) |>
    ggplot2::theme_minimal() |>
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) |>
    ggplot2::labs(
      title = "Inventory Completeness (Chao2) by Country",
      x = "Country",
      y = "Chao2 Completeness"
    )

  ggplot2::ggsave(output_path, p, width = 12, height = 6, dpi = 150)
}
