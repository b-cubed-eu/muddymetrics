library(testthat)

test_that("generate_site_reports handles missing template", {
  expect_error(
    generate_site_reports(
      data_dir = "/nonexistent",
      output_dir = "/nonexistent",
      template_path = "/nonexistent.Rmd"
    ),
    "Template not found"
  )
})


test_that("get_site_list returns empty list for missing directory", {
  result <- get_site_list("/nonexistent")
  expect_type(result, "list")
  expect_equal(length(result), 0)
})


test_that("generate_summary_report creates output", {
  skip_if_not(rmarkdown::pandoc_available(), "Pandoc not available")

  test_data <- data.frame(
    site_id = c("site1", "site2"),
    site_name = c("Site A", "Site B"),
    country = c("France", "Germany"),
    continent = c("Europe", "Europe"),
    density_km2 = c(0.5, 0.3),
    data_class = c("Data-Rich", "Data-Poor"),
    stringsAsFactors = FALSE
  )

  temp_file <- tempfile(fileext = ".html")

  expect_error(
    generate_summary_report(test_data, temp_file),
    NA
  )

  unlink(temp_file)
})
