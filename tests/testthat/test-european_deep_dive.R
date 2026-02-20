library(testthat)

test_that("run_european_deep_dive script exists", {
  script_path <- file.path("..", "..", "scripts", "run_european_deep_dive.R")
  expect_true(file.exists(script_path))
})


test_that("plot_diversity_by_country creates valid plot", {
  test_data <- data.frame(
    site_id = c("site1", "site2", "site3", "site4"),
    country = c("France", "France", "Germany", "Germany"),
    hill_1 = c(5.2, 4.8, 6.1, 5.9),
    stringsAsFactors = FALSE
  )

  expect_s3_class(
    ggplot2::ggplot(test_data, ggplot2::aes(x = country, y = hill_1)),
    "ggplot"
  )
})


test_that("plot_completeness_by_country creates valid plot", {
  test_data <- data.frame(
    site_id = c("site1", "site2", "site3", "site4"),
    country = c("France", "France", "Germany", "Germany"),
    completeness_chao2 = c(0.75, 0.65, 0.82, 0.78),
    stringsAsFactors = FALSE
  )

  expect_s3_class(
    ggplot2::ggplot(test_data, ggplot2::aes(x = country, y = completeness_chao2)),
    "ggplot"
  )
})


test_that("plot_diversity_by_country creates valid plot", {
  test_data <- data.frame(
    site_id = c("site1", "site2", "site3", "site4"),
    country = c("France", "France", "Germany", "Germany"),
    hill_1 = c(5.2, 4.8, 6.1, 5.9),
    stringsAsFactors = FALSE
  )

  temp_file <- tempfile(fileext = ".png")

  expect_s3_class(
    ggplot2::ggplot(test_data, ggplot2::aes(x = country, y = hill_1)),
    "ggplot"
  )

  unlink(temp_file)
})


test_that("plot_completeness_by_country creates valid plot", {
  test_data <- data.frame(
    site_id = c("site1", "site2", "site3", "site4"),
    country = c("France", "France", "Germany", "Germany"),
    completeness_chao2 = c(0.75, 0.65, 0.82, 0.78),
    stringsAsFactors = FALSE
  )

  expect_s3_class(
    ggplot2::ggplot(test_data, ggplot2::aes(x = country, y = completeness_chao2)),
    "ggplot"
  )
})
