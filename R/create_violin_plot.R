#' Create Violin Plot of Site Occurrences
#'
#' @param nested_list A nested list of countries and sites with occurrence values.
#' @return A ggplot object.
#' @importFrom purrr map map_dbl map_dfr set_names
#' @importFrom dplyr mutate filter group_by summarise arrange pull
#' @importFrom ggplot2 ggplot aes geom_violin geom_jitter stat_summary labs scale_y_log10 theme_minimal theme element_text element_blank
#' @importFrom tibble tibble
#' @importFrom rlang .data
#' @importFrom stats median runif rnorm
#' @importFrom utils head str
#' @export
create_violin_plot <- function(nested_list) {

  # 1. DATA TRANSFORMATION --------------------------------------------------

  # Get the list of country names from the input data
  country_names_from_list <- names(nested_list)

  # Flatten the nested list and clean values
  df_occurrences <- purrr::map_dfr(nested_list, function(.x) {
    site_names <- names(.x)
    raw_values <- unname(.x)

    # Process raw values: Handle NULL, strip units, and convert to numeric
    cleaned_values <- purrr::map_dbl(raw_values, function(val) {
      if (is.null(val) || (length(val) == 1 && is.na(val))) {
        return(NA_real_)
      }

      # Strip unit string and convert to numeric
      cleaned_str <- sub(" \\[1/km\\^2\\]", "", as.character(val), fixed = TRUE)
      return(as.numeric(cleaned_str))
    })

    # Create tibble for this country
    tibble::tibble(
      Value = cleaned_values,
      Site = site_names
    )
  }, .id = "Country") |>
    # Use .data$ to prevent R CMD check "undefined global variable" notes
    dplyr::mutate(Country = factor(.data$Country, levels = country_names_from_list)) |>
    dplyr::filter(!is.na(.data$Value))

  # 2. REORDERING -----------------------------------------------------------

  # Reorder countries by median value for better readability
  country_order <- df_occurrences |>
    dplyr::group_by(.data$Country) |>
    dplyr::summarise(MedianValue = stats::median(.data$Value, na.rm = TRUE)) |>
    dplyr::arrange(dplyr::desc(.data$MedianValue)) |>
    dplyr::pull(.data$Country)

  df_occurrences$Country <- factor(df_occurrences$Country, levels = country_order)

  # 3. VISUALIZATION --------------------------------------------------------

  occurrence_plot <- ggplot2::ggplot(df_occurrences, ggplot2::aes(x = .data$Country, y = .data$Value)) +
    ggplot2::geom_violin(
      ggplot2::aes(fill = .data$Country),
      alpha = 0.6,
      show.legend = FALSE,
      width = 1.2
    ) +
    ggplot2::geom_jitter(
      color = "black",
      size = 1.5,
      alpha = 0.7,
      width = 0.25
    ) +
    ggplot2::stat_summary(
      fun.data = function(x) {
        y_med <- stats::median(x, na.rm = TRUE)
        # Small buffer for the pointrange geom
        data.frame(y = y_med, ymin = y_med - 0.1, ymax = y_med + 0.1)
      },
      geom = "pointrange",
      linewidth = 1,
      color = "darkred"
    ) +
    ggplot2::labs(
      title = "Site Occurrence Density Distribution by Country (Units Removed)",
      subtitle = paste0(
        "Countries (n=", length(country_order),
        ") ordered by median occurrence density | Total Sites: ",
        nrow(df_occurrences)
      ),
      y = "Median Density (Occurrences per km^2)",
      x = NULL
    ) +
    ggplot2::scale_y_log10(
      labels = scales::label_number()
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 16),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "gray40")
    )

  return(occurrence_plot)
}
