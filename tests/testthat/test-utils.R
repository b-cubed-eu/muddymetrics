library(testthat)
library(muddymetrics)

test_that("check_gbif_credentials works", {
  # Mocking options
  op <- options(gbif_user = NULL, gbif_pwd = NULL, gbif_email = NULL)
  on.exit(options(op))
  
  expect_error(check_gbif_credentials(), "GBIF credentials are not set")
  
  options(gbif_user = "user", gbif_pwd = "pwd", gbif_email = "email")
  expect_true(check_gbif_credentials())
})
