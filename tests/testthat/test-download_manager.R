library(testthat)
library(muddymetrics)

test_that("save_download_key and load_download_key work", {
  tmp_file <- tempfile()
  key <- "test-key-123"
  
  save_download_key(key, tmp_file)
  expect_true(file.exists(tmp_file))
  
  loaded_key <- load_download_key(tmp_file)
  expect_equal(loaded_key, key)
  
  if (file.exists(tmp_file)) file.remove(tmp_file)
})

test_that("load_download_key returns NULL if file missing", {
  expect_null(load_download_key("non_existent_file"))
})

test_that("download_robust constructs correct command", {
  # We can't easily test system() calls without mocking or checking output
  # But we can check it runs (even if it fails on invalid URL)
  expect_message(download_robust("invalid_url", "dest"), "Starting robust download")
  if (file.exists("dest")) file.remove("dest")
})
