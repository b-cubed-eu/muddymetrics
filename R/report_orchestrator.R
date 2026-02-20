#' Generate Site Reports in Batch
#'
#' @description
#' This script orchestrates the generation of individual site reports
#' by rendering the RMarkdown template for each site in the dataset.
#'
#' @param data_dir Character. Path to the base directory containing site data.
#' @param output_dir Character. Path to the output directory with indicator results.
#' @param template_path Character. Path to the RMarkdown template.
#' @param reports_dir Character. Directory to save generated reports.
#' @param format Character. Output format: "html", "pdf", or "word".
#'
#' @export
generate_site_reports <- function(data_dir,
                                   output_dir,
                                   template_path = NULL,
                                   reports_dir = "output/reports",
                                   format = "html") {

  if (is.null(template_path)) {
    template_path <- system.file("templates", "site_report_template.Rmd",
                                 package = "muddymetrics")
  }

  if (!file.exists(template_path)) {
    stop("Template not found: ", template_path)
  }

  if (!dir.exists(reports_dir)) {
    dir.create(reports_dir, recursive = TRUE)
  }

  sites <- get_site_list(data_dir)

  message(paste0("Generating reports for ", length(sites), " sites..."))

  for (i in seq_along(sites)) {
    site_info <- sites[[i]]

    tryCatch({
      message(paste0("[", i, "/", length(sites), "] Generating: ", site_info$site_name))

      output_file <- file.path(
        reports_dir,
        paste0(site_info$site_id, "_report.", if (format == "html") "html" else format)
      )

      rmarkdown::render(
        input = template_path,
        output_file = output_file,
        params = list(
          site_id = site_info$site_id,
          site_name = site_info$site_name,
          country = site_info$country,
          data_dir = data_dir,
          output_dir = output_dir
        ),
        quiet = TRUE
      )

    }, error = function(e) {
      message(paste0("Error generating report for ", site_info$site_id, ": ",
                     conditionMessage(e)))
    })
  }

  message(paste0("\nReports generated in: ", reports_dir))
  message(paste0("Total: ", length(sites), " sites"))
}


#' Get Site List from Data Directory
#'
#' @param data_dir Character. Path to data directory.
#'
#' @return A list of site information.
#'
get_site_list <- function(data_dir) {

  continents <- c("africa", "antarctica", "asia", "europe", "northamerica", "oceania", "southamerica")

  sites <- list()

  for (continent in continents) {
    continent_dir <- file.path(data_dir, continent)

    if (!dir.exists(continent_dir)) next

    countries <- list.files(continent_dir)

    for (country in countries) {
      country_dir <- file.path(continent_dir, country)

      if (!dir.exists(country_dir)) next

      files <- list.files(country_dir, pattern = "\\.csv$")

      for (f in files) {
        site_name <- sub("_data\\.csv$", "", f)
        site_id <- paste0(country, "_", site_name)

        sites[[length(sites) + 1]] <- list(
          site_id = site_id,
          site_name = site_name,
          country = country,
          continent = continent
        )
      }
    }
  }

  return(sites)
}


#' Generate Summary Dashboard Report
#'
#' @param summary_data Data frame. The aggregated global summary data.
#' @param output_path Character. Path to save the report.
#'
#' @export
generate_summary_report <- function(summary_data, output_path) {

  report_template <- "
# Global Ramsar Data Sufficiency Summary

## Overview

- **Total Sites Analyzed**: `r nrow(summary_data)`
- **Data-Rich Sites**: `r sum(summary_data$data_class == 'Data-Rich', na.rm = TRUE)`
- **Data-Poor Sites**: `r sum(summary_data$data_class == 'Data-Poor', na.rm = TRUE)`

## Summary Statistics

| Metric | Mean | Median | Min | Max |
|--------|------|--------|-----|-----|
| Density (records/km²) | `r mean(summary_data$density_km2, na.rm = TRUE)` | `r median(summary_data$density_km2, na.rm = TRUE)` | `r min(summary_data$density_km2, na.rm = TRUE)` | `r max(summary_data$density_km2, na.rm = TRUE)` |

## By Continent

```{r continent-summary, echo=FALSE}
summary_data |>
  dplyr::group_by(continent) |>
  dplyr::summarise(
    n = dplyr::n(),
    data_rich = sum(data_class == 'Data-Rich', na.rm = TRUE),
    pct = round(data_rich / n * 100, 1)
  ) |>
  knitr::kable()
```

Generated: `r Sys.Date()`
"

  temp_file <- tempfile(fileext = ".Rmd")
  writeLines(report_template, temp_file)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  rmarkdown::render(
    input = temp_file,
    output_file = output_path,
    quiet = TRUE
  )

  unlink(temp_file)

  message(paste0("Summary report saved to: ", output_path))
}
