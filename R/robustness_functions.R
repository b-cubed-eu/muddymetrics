#' Assess Indicator Robustness Using Dubicube
#'
#' @description
#' Wrapper function to assess the robustness of biodiversity indicators using
#' the dubicube package. Uses bootstrap resampling to estimate confidence
#' intervals and uncertainty metrics.
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param indicator_function Function. A function that calculates the indicator.
#' @param grouping_var Character. The grouping variable for bootstrap (e.g., "year").
#' @param n_samples Numeric. Number of bootstrap samples (default 100).
#' @param conf_level Numeric. Confidence level for intervals (default 0.95).
#' @param method Character. Method for CI calculation: "perc", "bca", "norm", "basic".
#' @param shapefile_path Optional path to WKT file for spatial filtering.
#'
#' @return A list with robustness metrics including:
#'   - original_estimate: The original indicator value
#'   - bootstrap_se: Standard error from bootstrap
#'   - ci_lower: Lower confidence interval bound
#'   - ci_upper: Upper confidence interval bound
#'   - bias: Estimated bias
#'   - cv_result: Cross-validation results (if available)
#'
#' @export
assess_indicator_robustness <- function(cube,
                                         indicator_function,
                                         grouping_var = "year",
                                         n_samples = 100,
                                         conf_level = 0.95,
                                         method = "perc",
                                         shapefile_path = NULL) {

  if (!requireNamespace("dubicube", quietly = TRUE)) {
    stop("dubicube package required. Install with: install.packages('dubicube', repos = 'https://b-cubed-eu.r-universe.dev')")
  }

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  cube_processed <- tryCatch({
    b3gbi::process_cube(cube, separator = ",")
  }, error = function(e) {
    message("Error processing cube: ", conditionMessage(e))
    NULL
  })

  if (is.null(cube_processed)) {
    cube_processed <- list(data = cube)
  }

  original_estimate <- tryCatch({
    indicator_function(cube_processed$data)
  }, error = function(e) {
    message("Error calculating indicator: ", conditionMessage(e))
    NA
  })

  bootstrap_result <- tryCatch({
    dubicube::bootstrap_cube(
      data_cube = cube_processed,
      fun = indicator_function,
      grouping_var = grouping_var,
      samples = n_samples,
      processed_cube = TRUE
    )
  }, error = function(e) {
    message("Bootstrap error: ", conditionMessage(e))
    NULL
  })

  if (is.null(bootstrap_result)) {
    return(list(
      original_estimate = original_estimate,
      bootstrap_se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      bias = NA_real_,
      method = method,
      n_samples = n_samples
    ))
  }

  bootstrap_df <- dubicube::boot_list_to_dataframe(bootstrap_result)

  ci_result <- tryCatch({
    dubicube::calculate_bootstrap_ci(
      bootstrap_samples_df = bootstrap_df,
      grouping_var = grouping_var,
      type = method,
      conf = conf_level
    )
  }, error = function(e) {
    message("CI calculation error: ", conditionMessage(e))
    NULL
  })

  if (is.null(ci_result)) {
    return(list(
      original_estimate = original_estimate,
      bootstrap_se = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      bias = NA_real_,
      method = method,
      n_samples = n_samples
    ))
  }

  t_star <- bootstrap_df$diversity_val
  t0 <- original_estimate
  bias_est <- mean(t_star) - t0
  bootstrap_se <- stats::sd(t_star)

  return(list(
    original_estimate = original_estimate,
    bootstrap_se = bootstrap_se,
    ci_lower = ci_result$ci_lower,
    ci_upper = ci_result$ci_upper,
    bias = bias_est,
    method = method,
    conf_level = conf_level,
    n_samples = n_samples,
    bootstrap_distribution = t_star
  ))
}


#' Cross-Validate Indicator
#'
#' @description
#' Performs leave-one-out cross-validation on an indicator to assess its stability.
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param indicator_function Function. A function that calculates the indicator.
#' @param grouping_var Character. The grouping variable for CV (e.g., "year").
#' @param shapefile_path Optional path to WKT file for spatial filtering.
#'
#' @return A list with CV metrics.
#'
#' @export
cross_validate_indicator <- function(cube,
                                     indicator_function,
                                     grouping_var = "year",
                                     shapefile_path = NULL) {

  if (!requireNamespace("dubicube", quietly = TRUE)) {
    stop("dubicube package required.")
  }

  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }

  if (!is.null(shapefile_path) && file.exists(shapefile_path)) {
    cube <- filter_by_shapefile(cube, shapefile_path)
  }

  cube_processed <- tryCatch({
    b3gbi::process_cube(cube, separator = ",")
  }, error = function(e) {
    list(data = cube)
  })

  cv_result <- tryCatch({
    dubicube::cross_validate_cube(
      data_cube = cube_processed,
      fun = indicator_function,
      grouping_var = grouping_var,
      processed_cube = TRUE
    )
  }, error = function(e) {
    message("CV error: ", conditionMessage(e))
    NULL
  })

  if (is.null(cv_result)) {
    return(list(
      cv_rmse = NA_real_,
      cv_mae = NA_real_,
      n_removed = 0
    ))
  }

  cv_rmse <- sqrt(mean(cv_result$sq_error, na.rm = TRUE))
  cv_mae <- mean(cv_result$abs_error, na.rm = TRUE)

  return(list(
    cv_rmse = cv_rmse,
    cv_mae = cv_mae,
    n_removed = sum(cv_result$removed),
    cv_details = cv_result
  ))
}


#' Batch Assess Robustness for Multiple Sites
#'
#' @param inputdir Character. Path to directory containing site CSV files.
#' @param outputdir Character. Path to save results.
#' @param shapefiledir Character. Path to WKT files.
#' @param indicator_function Function. The indicator function to assess.
#' @param grouping_var Character. Grouping variable for bootstrap.
#' @param n_samples Numeric. Number of bootstrap samples.
#'
#' @return A data frame with robustness metrics per site.
#'
#' @export
assess_robustness_batch <- function(inputdir,
                                    outputdir,
                                    shapefiledir,
                                    indicator_function,
                                    grouping_var = "year",
                                    n_samples = 100) {

  if (!dir.exists(outputdir)) {
    dir.create(outputdir, recursive = TRUE)
  }

  results <- data.frame(
    site_id = character(),
    country = character(),
    original_estimate = numeric(),
    bootstrap_se = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    bias = numeric(),
    cv_rmse = numeric(),
    cv_mae = numeric(),
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

        robustness <- assess_indicator_robustness(
          cube = cubepath,
          indicator_function = indicator_function,
          grouping_var = grouping_var,
          n_samples = n_samples,
          shapefile_path = shapefilepath
        )

        cv_result <- cross_validate_indicator(
          cube = cubepath,
          indicator_function = indicator_function,
          grouping_var = grouping_var,
          shapefile_path = shapefilepath
        )

        results <- rbind(results, data.frame(
          site_id = site_id,
          country = country_name,
          original_estimate = robustness$original_estimate,
          bootstrap_se = robustness$bootstrap_se,
          ci_lower = robustness$ci_lower,
          ci_upper = robustness$ci_upper,
          bias = robustness$bias,
          cv_rmse = cv_result$cv_rmse,
          cv_mae = cv_result$cv_mae,
          stringsAsFactors = FALSE
        ))

      }, error = function(e) {
        message(paste0("Error processing ", site_file, ": ", conditionMessage(e)))
      })
    }
  }

  return(results)
}
