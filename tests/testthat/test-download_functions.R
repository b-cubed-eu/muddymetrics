library(testthat)

test_that("get_gbif_predicates returns correct predicates", {
  # Mocking gbif_base_filters
  filters <- list(
    hasGeospatialIssue = FALSE,
    hasCoordinate = TRUE,
    basisOfRecord = c("OBSERVATION", "HUMAN_OBSERVATION")
  )
  
  wkt <- "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))"
  
  # This function doesn't exist yet, which is the point of TDD (Red phase)
  preds <- get_gbif_predicates(wkt, filters)
  
  expect_s3_class(preds, "occ_predicate_list")
})
