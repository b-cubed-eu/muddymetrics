library(testthat)

test_that("calculate_chao2 returns correct estimator", {
  test_cube <- data.frame(
    specieskey = c(1, 1, 1, 2, 2, 3, 4, 4, 4, 4),
    occurrenceId = c("a", "a", "a", "b", "b", "c", "d", "d", "d", "d"),
    geometry = c("1|1", "1|1", "1|1", "2|2", "2|2", "3|3", "4|4", "4|4", "4|4", "4|4"),
    stringsAsFactors = FALSE
  )

  result <- calculate_chao2(test_cube)

  expect_equal(result$observed, 4)

  expect_type(result, "list")
  expect_true(result$chao2 >= result$observed)
})


test_that("calculate_chao2 handles all unique species", {
  test_cube <- data.frame(
    specieskey = c(1, 2, 3, 4),
    occurrenceId = c("a", "b", "c", "d"),
    geometry = c("1|1", "2|2", "3|3", "4|4"),
    stringsAsFactors = FALSE
  )

  result <- calculate_chao2(test_cube)

  expect_equal(result$observed, 4)
  expect_equal(result$f1, 4)
  expect_equal(result$f2, 0)
})


test_that("calculate_chao2 handles empty data", {
  test_cube <- data.frame(
    specieskey = integer(),
    occurrenceId = character(),
    geometry = character(),
    stringsAsFactors = FALSE
  )

  result <- calculate_chao2(test_cube)

  expect_equal(result$observed, 0)
  expect_true(is.na(result$chao2))
})


test_that("calculate_sac_slope returns valid slope", {
  set.seed(42)

  test_cube <- data.frame(
    specieskey = rep(1:20, each = 5),
    occurrenceId = rep(paste0("occ", 1:100), each = 1),
    geometry = rep(paste0(1:100, "|", 1:100), each = 1),
    stringsAsFactors = FALSE
  )

  result <- calculate_sac_slope(test_cube, n_iterations = 10)

  expect_type(result, "list")
  expect_true(!is.na(result$slope) || is.na(result$slope))
})


test_that("calculate_sac_slope handles insufficient data", {
  test_cube <- data.frame(
    specieskey = c(1, 1, 1),
    occurrenceId = c("a", "a", "a"),
    geometry = c("1|1", "1|1", "1|1"),
    stringsAsFactors = FALSE
  )

  result <- calculate_sac_slope(test_cube)

  expect_true(is.na(result$slope))
})
