library(testthat)
library(b3gbi)
library(ggplot2)

test_that("save_ramsar_plot saves a ggplot object", {
  # Create a dummy indicator object
  data("example_cube_1", package = "b3gbi")
  # Source needed file for helper function if needed, but here we can just mock or use the real one if sourced
  # For the test, we need calculate_ramsar_metric
  source("../../R/indicator_functions.R")
  
  ind <- calculate_ramsar_metric(example_cube_1, "obs_richness", "ts", ci_type = "none")
  
  temp_file <- tempfile(fileext = ".png")
  
  # This function doesn't exist yet (Red phase)
  save_ramsar_plot(ind, temp_file)
  
  expect_true(file.exists(temp_file))
  if (file.exists(temp_file)) file.remove(temp_file)
})
