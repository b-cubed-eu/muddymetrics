library(testthat)

test_that("calculate_mean_year returns correct weighted mean", {
  test_cube <- data.frame(
    year = c(2020, 2019, 2018, 2021),
    occurrences = c(100, 50, 25, 75),
    geometry = c("1|1", "2|2", "3|3", "4|4"),
    stringsAsFactors = FALSE
  )

  result <- calculate_mean_year(test_cube)

  expected_mean <- (2020 * 100 + 2019 * 50 + 2018 * 25 + 2021 * 75) / (100 + 50 + 25 + 75)

  expect_equal(result, expected_mean)
})


test_that("calculate_mean_year handles NA values", {
  test_cube <- data.frame(
    year = c(2020, NA, 2019, 2021),
    occurrences = c(100, 50, 25, 75),
    geometry = c("1|1", "2|2", "3|3", "4|4"),
    stringsAsFactors = FALSE
  )

  result <- calculate_mean_year(test_cube)

  expected_mean <- (2020 * 100 + 2019 * 25 + 2021 * 75) / (100 + 25 + 75)

  expect_equal(result, expected_mean)
})


test_that("calculate_mean_year returns NA for all invalid data", {
  test_cube <- data.frame(
    year = c(NA, NA, NA),
    occurrences = c(NA, NA, NA),
    geometry = c("1|1", "2|2", "3|3"),
    stringsAsFactors = FALSE
  )

  result <- calculate_mean_year(test_cube)

  expect_true(is.na(result))
})
