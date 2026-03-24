library(testthat)
library(b3gbi)

test_that("calc_ramsar_indicator processes a real site flow", {
  # Locate the internal test data
  test_base <- system.file("extdata/test_data", package = "muddymetrics")

  # Fallback for local development if system.file fails
  if (test_base == "") test_base <- "inst/extdata/test_data"

  skip_if(!dir.exists(test_base), "Test data directory not found")

  # Setup temp output
  tmp_main <- file.path(tempdir(), "ramsar_output")
  if (dir.exists(tmp_main)) unlink(tmp_main, recursive = TRUE)
  dir.create(tmp_main)

  # Run the actual workflow
  # Using a try-catch in the test to capture b3gbi internal failures
  expect_error(
    result <- calc_ramsar_indicator(
      indicator = "obs_richness_ts",
      inputdir = file.path(test_base, "input"),
      maindir = tmp_main,
      shapefiledir = file.path(test_base, "wkt"),
      continent = "Asia",
      plot_args = list(smoothed_trend = FALSE)
    ),
    regexp = NA # Expect no error
  )

  # Validate results
  expect_type(result, "list")
  expect_true("Vietnam" %in% names(result$mean))

  # Check if ggsave actually produced the PNG
  out_files <- list.files(tmp_main, recursive = TRUE, pattern = "\\.png$")
  expect_gt(length(out_files), 0)

  # Cleanup
  unlink(tmp_main, recursive = TRUE)
})
