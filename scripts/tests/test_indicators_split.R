# Tests for the split interval runs (run_indicators.R --split-hours): the parts,
# merged, must give exactly the results of one unsplit run.
# Run from the repository root:
#   Rscript -e 'testthat::test_file("scripts/tests/test_indicators_split.R")'

library(testthat)
suppressPackageStartupMessages({ library(data.table) })
here <- if (file.exists("scripts/lib/indicators_v1.R")) "scripts" else ".."
source(file.path(here, "lib/indicators_v1.R"))

test_that("the job plan splits by the earlier run's time, longest parts first", {
  prior <- data.table(site_id = c("s1", "s1", "s2", "s2"), subset = "aves", grid = "1km", level = c(0.9, 0.95, 0.9, 0.95),
                      seconds = c(3600, 3600, 60, 60), year = 2000, q = 0)
  p <- ci_job_plan(data.table(site_id = c("s1", "s2", "s3"), subset = "aves"), prior = prior, hours = 0.5, default_k = 4)
  expect_equal(p[site_id == "s1", unique(k)], 6L)         # 2 h x 1.3 / 0.5 h
  expect_equal(p[site_id == "s2", unique(k)], 1L)
  expect_equal(p[site_id == "s3", unique(k)], 4L)         # no earlier time
  expect_equal(p[site_id == "s1", sort(years)], sprintf("m6r%d", 0:5))
  expect_equal(p$site_id[1], "s3")                        # 0.5 h parts (default) before 0.43 h parts
})

test_that("combinations without an earlier time get a time predicted from their size", {
  set.seed(2)
  n <- 30
  sz <- data.table(site_id = sprintf("s%d", 1:(n + 1)), subset = "aves", n_rows = round(10^runif(n + 1, 3, 6)),
                   n_cells = round(10^runif(n + 1, 1, 3)), S_obs = round(10^runif(n + 1, 1.5, 3)))
  sz[, h := 1e-4 * n_rows^0.5 * n_cells^0.4 * S_obs^0.3]
  prior <- sz[1:n, .(site_id, subset, grid = "1km", level = 0.9, seconds = h * 3600, year = 2000, q = 0)]
  sz[n + 1, `:=`(n_rows = 1e6, n_cells = 1000, S_obs = 1000)]
  p <- ci_job_plan(sz[, .(site_id, subset)], prior = prior, hours = 0.5, size = sz, growth = 1)
  u <- unique(p[, .(site_id, k, est_hours, predicted)])
  new <- u[site_id == sprintf("s%d", n + 1)]
  expect_true(new$predicted)
  expect_equal(new$est_hours, 1.5 * 1e-4 * 1e6^0.5 * 1000^0.4 * 1000^0.3, tolerance = 0.01)
  expect_false(any(u$predicted[u$site_id != new$site_id]))
})

skip_if_not_installed("b3gbi")
skip_if_not_installed("withr")

test_that("a split interval run gives the same results as an unsplit one", {
  repo <- normalizePath(file.path(here, ".."))
  d <- tempfile("split_"); dir.create(d)
  file.symlink(file.path(repo, "scripts"), file.path(d, "scripts"))
  cdir <- file.path(d, "inst/extdata/ramsar_site_data_100m_europe/Belgium"); dir.create(cdir, recursive = TRUE)
  set.seed(1)
  mk <- function(n_sp, n_cells, years) {
    g <- CJ(year = years, cell = sprintf("31UES%02d0%02d0", 1:n_cells, 1:n_cells)[seq_len(n_cells)], sp = seq_len(n_sp))
    g <- g[runif(.N) < 0.25]
    data.table(kingdom = "Animalia", kingdomkey = "1", phylum = "Chordata", phylumkey = "44", class = "Aves", classkey = "212",
               order = "O", orderkey = "1", family = "F", familykey = "1", genus = "G", genuskey = "1",
               species = paste("Avis", g$sp), specieskey = as.character(1000 + g$sp), year = g$year, mgrscellcode = g$cell,
               kingdomcount = 1, phylumcount = 1, classcount = 1, ordercount = 1, familycount = 1, genuscount = 1,
               occurrences = 1, mintemporaluncertainty = 60, mincoordinateuncertaintyinmeters = 10)
  }
  fwrite(mk(25, 30, 2015:2022), file.path(cdir, "site_1_One_data.csv"))
  fwrite(mk(15, 20, 2018:2021), file.path(cdir, "site_2_Two_data.csv"))
  # a screen file in which both sites pass for birds
  # (other sites without cubes, all failing, so that the temporal elbow is defined)
  ids <- c("site_1", "site_2", sprintf("site_%d", 3:20))
  raw <- data.table(continent = "europe", country = "Belgium", site_id = ids, subset = "aves", grid = "1km",
                    status = "computed", density_km2 = c(10, 10, rep(0, 18)), chao2_completeness = 0.9,
                    temporal_cha = c(1, 1, seq(0.05, 0.9, length.out = 18)),
                    precision_ratio = 0.01, precision_ratio_median = 0.01, precision_ratio_d61 = 0.01,
                    coverage_deficit_min_cells = 0, coverage_deficit_b3gbi = 0, coverage_deficit_informative = 0,
                    n_years_informative = 20, mean_cov_min_cells = 1, mean_cov_b3gbi = 1, mean_cov_informative = 1)
  fwrite(raw, file.path(d, "raw.csv"))
  run <- function(...) withr::with_dir(d, system2("Rscript", c("scripts/run_indicators.R", "--screen", "raw.csv", "--mode", "ci",
                                                               "--nboot", "20", "--grids", "1km", "--cores", "2", ...),
                                                  stdout = TRUE, stderr = TRUE))
  a <- run("--out", "out_a"); expect_null(attr(a, "status"))
  b <- run("--out", "out_b", "--split-hours", "0.5", "--split-default", "3"); expect_null(attr(b, "status"))
  fa <- fread(file.path(d, "out_a/hill_ci_1km_b20.csv")); fb <- fread(file.path(d, "out_b/hill_ci_1km_b20.csv"))
  expect_gt(nrow(fa), 20)
  expect_length(list.files(file.path(d, "out_b/parts"), pattern = "_ym3r[0-2]\\.csv$"), 6)
  expect_equal(names(fa), names(fb))
  cols <- setdiff(names(fa), "seconds")
  setorderv(fa, c("site_id", "subset", "grid", "level", "year", "q")); setorderv(fb, c("site_id", "subset", "grid", "level", "year", "q"))
  expect_equal(fa[, ..cols], fb[, ..cols])
  expect_true(any(fa$ci_status == "ok"))
})
