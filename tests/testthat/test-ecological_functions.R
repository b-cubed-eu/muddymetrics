library(testthat)

test_that("calculate_alpha_diversity returns correct Hill numbers", {
  test_cube <- data.frame(
    specieskey = c(1, 1, 1, 2, 2, 3, 4, 4, 4, 4),
    occurrences = c(100, 100, 100, 50, 50, 25, 10, 10, 10, 10),
    geometry = rep("1|1", 10),
    stringsAsFactors = FALSE
  )

  result <- calculate_alpha_diversity(test_cube)

  expect_equal(result$hill_0, 4)
  expect_type(result, "list")
  expect_true(result$hill_1 <= result$hill_0 + 1)
  expect_true(result$hill_2 <= result$hill_1)
})


test_that("calculate_alpha_diversity handles equal abundances", {
  test_cube <- data.frame(
    specieskey = c(1, 2, 3, 4),
    occurrences = c(25, 25, 25, 25),
    geometry = rep("1|1", 4),
    stringsAsFactors = FALSE
  )

  result <- calculate_alpha_diversity(test_cube)

  expect_equal(result$hill_0, 4)
  expect_equal(result$hill_1, 4, tolerance = 0.01)
})


test_that("calculate_alpha_diversity handles empty data", {
  test_cube <- data.frame(
    specieskey = integer(),
    occurrences = integer(),
    geometry = character(),
    stringsAsFactors = FALSE
  )

  result <- calculate_alpha_diversity(test_cube)

  expect_equal(result$hill_0, NA_integer_)
})


test_that("calculate_evenness returns valid Pielou J", {
  test_cube <- data.frame(
    specieskey = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3),
    occurrences = c(50, 50, 50, 30, 30, 30, 10, 10, 10, 10),
    geometry = rep("1|1", 10),
    stringsAsFactors = FALSE
  )

  result <- calculate_evenness(test_cube)

  expect_type(result, "list")
  expect_true(result$pielou_j >= 0 && result$pielou_j <= 1)
})


test_that("calculate_evenness handles single species", {
  test_cube <- data.frame(
    specieskey = c(1, 1, 1),
    occurrences = c(10, 10, 10),
    geometry = rep("1|1", 3),
    stringsAsFactors = FALSE
  )

  result <- calculate_evenness(test_cube)

  expect_equal(result$pielou_j, NA_real_)
})


test_that("calculate_rarity_metrics returns valid metrics", {
  test_cube <- data.frame(
    specieskey = c(1, 1, 1, 2, 2, 3, 3, 3, 3, 3),
    occurrences = c(100, 100, 100, 10, 10, 5, 5, 5, 5, 5),
    geometry = c(rep("1|1", 3), rep("2|2", 2), rep("3|3", 5)),
    stringsAsFactors = FALSE
  )

  result <- calculate_rarity_metrics(test_cube, site_area_km2 = 100)

  expect_type(result, "list")
  expect_true(result$abundance_rarity >= 0 && result$abundance_rarity <= 1)
  expect_true(result$area_rarity >= 0)
})


test_that("calculate_rarity_metrics handles edge cases", {
  test_cube <- data.frame(
    specieskey = integer(),
    occurrences = integer(),
    geometry = character(),
    stringsAsFactors = FALSE
  )

  result <- calculate_rarity_metrics(test_cube)

  expect_equal(result$area_rarity, NA_real_)
})
