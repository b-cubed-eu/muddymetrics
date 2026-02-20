library(testthat)

test_that("generate_global_summary classifies sites correctly", {
  test_data <- data.frame(
    site_id = c("site1", "site2", "site3", "site4"),
    site_name = c("Site A", "Site B", "Site C", "Site D"),
    country = c("CountryA", "CountryA", "CountryB", "CountryB"),
    continent = c("Europe", "Europe", "Asia", "Asia"),
    obs_richness = c(100, 50, 75, 25),
    total_occurrences = c(1000, 100, 500, 10),
    cum_richness = c(120, 60, 80, 30),
    density_km2 = c(0.5, 0.1, 0.3, 0.05),
    mean_year = c(2020, 2015, 2018, 2010),
    year_min = c(2000, 2000, 2005, 2000),
    year_max = c(2023, 2022, 2023, 2020),
    chao2 = c(150, 80, 100, 40),
    chao2_completeness = c(0.80, 0.62, 0.75, 0.50),
    sac_slope = c(0.05, 0.15, 0.08, 0.20),
    stringsAsFactors = FALSE
  )

  temp_file <- tempfile(fileext = ".csv")

  result <- generate_global_summary(
    test_data,
    temp_file,
    threshold_density = 0.25,
    threshold_chao2 = 0.70,
    threshold_slope = 0.10
  )

  expect_true("data_class" %in% names(result))
  expect_equal(result$data_class[1], "Data-Rich")
  expect_equal(result$data_class[2], "Data-Poor")
  expect_true(result$passes_density[1])
  expect_false(result$passes_density[2])

  unlink(temp_file)
})


test_that("perform_data_gap_analysis returns correct statistics", {
  test_data <- data.frame(
    site_id = c("site1", "site2", "site3", "site4"),
    site_name = c("Site A", "Site B", "Site C", "Site D"),
    country = c("CountryA", "CountryA", "CountryB", "CountryB"),
    continent = c("Europe", "Europe", "Asia", "Asia"),
    obs_richness = NA,
    total_occurrences = NA,
    cum_richness = NA,
    density_km2 = c(0.5, 0.1, 0.3, 0.05),
    mean_year = NA,
    year_min = NA,
    year_max = NA,
    chao2 = NA,
    chao2_completeness = c(0.80, 0.62, 0.75, 0.50),
    sac_slope = c(0.05, 0.15, 0.08, 0.20),
    data_class = c("Data-Rich", "Data-Poor", "Data-Rich", "Data-Poor"),
    passes_density = c(TRUE, FALSE, TRUE, FALSE),
    passes_chao2 = c(TRUE, FALSE, TRUE, FALSE),
    passes_slope = c(TRUE, FALSE, TRUE, FALSE),
    passes_all_thresholds = c(TRUE, FALSE, FALSE, FALSE),
    troia_moderate = c(TRUE, FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  result <- perform_data_gap_analysis(test_data)

  expect_equal(result$total_sites, 4)
  expect_equal(result$data_rich_sites, 2)
  expect_equal(result$data_poor_sites, 2)
  expect_equal(result$pass_rate, 0.5)

  expect_true("by_continent" %in% names(result))
  expect_true("by_country" %in% names(result))
})
