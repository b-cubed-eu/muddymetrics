library(testthat)

test_that("calculate_invasive_impact handles missing EICAT data", {
  test_cube <- data.frame(
    specieskey = c(1, 2, 3),
    species = c("Species A", "Species B", "Species C"),
    occurrences = c(10, 20, 30),
    geometry = c("1|1", "2|2", "3|3"),
    stringsAsFactors = FALSE
  )

  eicat_empty <- data.frame(
    species = character(),
    impact_category = character(),
    mechanism = character(),
    stringsAsFactors = FALSE
  )

  result <- calculate_invasive_impact(test_cube, eicat_empty)

  expect_type(result, "list")
  expect_equal(result$total_invasive_species, 0)
})


test_that("calculate_invasive_impact with valid EICAT data", {
  test_cube <- data.frame(
    specieskey = c(1, 2, 3),
    species = c("Acacia melanoxylon", "Species B", "Species C"),
    occurrences = c(10, 20, 30),
    geometry = c("1|1", "2|2", "3|3"),
    stringsAsFactors = FALSE
  )

  eicat_data <- data.frame(
    species = c("Acacia melanoxylon"),
    impact_category = c("Moderate"),
    mechanism = c("Competition"),
    stringsAsFactors = FALSE
  )

  result <- calculate_invasive_impact(test_cube, eicat_data)

  expect_type(result, "list")
  expect_equal(result$total_species, 3)
})


test_that("calculate_phylogenetic_diversity checks for package", {
  test_cube <- data.frame(
    specieskey = c(1, 2, 3),
    occurrences = c(10, 20, 30),
    geometry = c("1|1", "2|2", "3|3"),
    stringsAsFactors = FALSE
  )

  expect_error(
    calculate_phylogenetic_diversity(test_cube, tree_path = "nonexistent.tre"),
    "not found"
  )
})


test_that("calculate_specialized_batch handles missing files gracefully", {
  results <- calculate_specialized_batch(
    inputdir = "/nonexistent",
    outputdir = tempdir(),
    shapefiledir = "/nonexistent"
  )

  expect_type(results, "list")
  expect_equal(nrow(results), 0)
})
