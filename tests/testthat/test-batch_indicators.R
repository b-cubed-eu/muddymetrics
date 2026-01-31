library(testthat)

test_that("calc_ramsar_indicator works with Vietnam test data", {
  # We'll use the existing Vietnam data for a small integration test
  # This avoids complex mocking of b3gbi internal behavior
  
  inputdir <- "inst/extdata/ramsar_site_data_100m_africa" # We'll use a subset or a fake one
  # Actually, let's create a temporary structure that mimics the Vietnam one
  
  tmp_root <- tempdir()
  tmp_input <- file.path(tmp_root, "input")
  tmp_output <- file.path(tmp_root, "output")
  tmp_wkt <- file.path(tmp_root, "wkt")
  
  dir.create(file.path(tmp_input, "Vietnam"), recursive = TRUE)
  dir.create(file.path(tmp_wkt, "Vietnam"), recursive = TRUE)
  dir.create(tmp_output)
  
  print(paste("Current WD:", getwd()))
  
  # Copy Vietnam data
  src_csv <- "../../inst/extdata/ramsar_site_data_100m_asia/Vietnam/site_2227_Lang_Sen_Wetland_Reserve_data.csv"
  dest_csv <- file.path(tmp_input, "Vietnam", "site_2227_Lang_Sen_Wetland_Reserve_data.csv")
  res_csv <- file.copy(src_csv, dest_csv)
  
  src_wkt <- "../../inst/extdata/ramsar_sites_wkt/Vietnam/site_2227_Lang_Sen_Wetland_Reserve.wkt"
  dest_wkt <- file.path(tmp_wkt, "Vietnam", "site_2227_Lang_Sen_Wetland_Reserve.wkt")
  res_wkt <- file.copy(src_wkt, dest_wkt)
  
  result <- calc_ramsar_indicator(
    indicator = "obs_richness_ts",
    inputdir = tmp_input,
    maindir = tmp_output,
    shapefiledir = tmp_wkt,
    continent = "Asia",
    plot_args = list(smoothed_trend = FALSE)
  )
  
  expect_type(result, "list")
  expect_true("Vietnam" %in% names(result$mean))
  expect_true("site_2227_Lang_Sen_Wetland_Reserve_data.csv" %in% names(result$mean$Vietnam))
  
  # Check if plot was saved
  expect_true(file.exists(file.path(tmp_output, "Vietnam", "site_2227_Lang_Sen_Wetland_Reserve_obs_richness_ts.png")))
  
  # Cleanup
  unlink(tmp_root, recursive = TRUE)
})