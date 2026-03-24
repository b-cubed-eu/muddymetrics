library(testthat)

test_that("assess_indicator_robustness handles missing package gracefully", {
  test_cube <- data.frame(
    year = rep(2020:2022, each = 10),
    specieskey = rep(1:10, 3),
    occurrences = runif(30, 1, 100),
    cellCode = rep(paste0("cell", 1:10), 3),
    stringsAsFactors = FALSE
  )

  simple_func <- function(data) {
    mean(data$occurrences, na.rm = TRUE)
  }

  skip_if_not_installed("dubicube")

  result <- assess_indicator_robustness(
    test_cube,
    indicator_function = simple_func,
    grouping_var = "year",
    n_samples = 10
  )

  expect_type(result, "list")
})


test_that("cross_validate_indicator handles errors gracefully", {
  test_cube <- data.frame(
    year = rep(2020:2022, each = 10),
    specieskey = rep(1:10, 3),
    occurrences = runif(30, 1, 100),
    cellCode = rep(paste0("cell", 1:10), 3),
    stringsAsFactors = FALSE
  )

  simple_func <- function(data) {
    mean(data$occurrences, na.rm = TRUE)
  }

  skip_if_not_installed("dubicube")

  result <- cross_validate_indicator(
    test_cube,
    indicator_function = simple_func,
    grouping_var = "year"
  )

  expect_type(result, "list")
})


test_that("assess_robustness_batch handles missing directories", {
  simple_func <- function(data) {
    mean(data$occurrences, na.rm = TRUE)
  }

  result <- assess_robustness_batch(
    inputdir = "/nonexistent",
    outputdir = tempdir(),
    shapefiledir = "/nonexistent",
    indicator_function = simple_func
  )

  expect_type(result, "list")
  expect_equal(nrow(result), 0)
})
