library(testthat)

test_that("validate_density_results returns correct summary statistics", {
  test_results <- data.frame(
    site_id = c("site1", "site2", "site3", "site4", "site5"),
    country = c("CountryA", "CountryA", "CountryB", "CountryB", "CountryC"),
    density_records_km2 = c(0.5, 0.2, 0.3, 0.1, 0.4),
    passes_threshold = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    stringsAsFactors = FALSE
  )

  result <- validate_density_results(test_results, threshold = 0.25)

  expect_type(result, "list")
  expect_equal(result$total_sites, 5)
  expect_equal(result$passing_sites, 3)
  expect_equal(result$failing_sites, 2)
  expect_equal(result$threshold, 0.25)
  expect_equal(result$pass_rate, 0.6)
  expect_equal(result$mean_density, 0.3)
  expect_equal(result$median_density, 0.3)
})


test_that("validate_density_results handles edge cases", {
  empty_results <- data.frame(
    site_id = character(),
    country = character(),
    density_records_km2 = numeric(),
    passes_threshold = logical(),
    stringsAsFactors = FALSE
  )

  result <- validate_density_results(empty_results, threshold = 0.25)

  expect_equal(result$total_sites, 0)
  expect_equal(result$passing_sites, 0)
  expect_true(is.nan(result$pass_rate))
})


test_that("validate_density_results with custom threshold", {
  test_results <- data.frame(
    site_id = c("site1", "site2"),
    country = c("CountryA", "CountryA"),
    density_records_km2 = c(0.3, 0.1),
    passes_threshold = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  result <- validate_density_results(test_results, threshold = 0.2)

  expect_equal(result$passing_sites, 1)
  expect_equal(result$threshold, 0.2)
})
