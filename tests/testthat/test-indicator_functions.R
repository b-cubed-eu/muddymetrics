library(testthat)
library(b3gbi)

test_that("calculate_ramsar_metric returns a valid indicator object", {
  # Use example data from b3gbi
  data("example_cube_1", package = "b3gbi")
  
  # This function doesn't exist yet (Red phase)
  result <- calculate_ramsar_metric(
    cube = example_cube_1,
    metric = "obs_richness",
    type = "ts",
    region = "Europe",
    ci_type = "none"
  )
  
  expect_s3_class(result, "indicator_ts")
  # b3gbi objects often store metadata in an 'attributes' or similar. 
  # Let's just check it has data.
  expect_true(!is.null(result$data))
})
