# Unit tests for scripts/lib/sufficiency_criteria.R
# Run from the repository root:
#   Rscript -e 'testthat::test_file("scripts/tests/test_sufficiency_criteria.R")'
#
# Expected values are either worked by hand (shown in comments) or checked
# against an independent reference implementation (vegan, iNEXT, b3gbi).

library(testthat)
# test_file() runs from this file's directory; also allow running from the repo root.
lib <- c("scripts/lib/sufficiency_criteria.R", "../lib/sufficiency_criteria.R")
source(lib[file.exists(lib)][1])

# A small cube: species x 100 m cell x year rows. Cells are valid MGRS codes.
make_cube <- function(rows) {
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(df) <- c("specieskey", "mgrscellcode", "year", "occurrences")
  df$year <- as.integer(df$year)
  df$occurrences <- as.numeric(df$occurrences)
  df$species <- paste0("sp", df$specieskey)
  df$kingdom <- "Animalia"; df$class <- "Aves"
  df$mincoordinateuncertaintyinmeters <- 50
  df
}
C <- c("32UNU001001", "32UNU002001", "32UNU003001", "32UNU004001", "32UNU005001")

# Occupancy (cells per species, years pooled): s1:1, s2:1, s3:2, s4:2, s5:3, s6:5
# -> S = 6, T = 5, Q1 = 2, Q2 = 2, U = 1+1+2+2+3+5 = 14
cube1 <- make_cube(list(
  c(1, C[1], 2000, 1),
  c(2, C[2], 2000, 1),
  c(3, C[1], 2000, 2), c(3, C[2], 2001, 1),
  c(4, C[3], 2001, 1), c(4, C[4], 2002, 4),
  c(5, C[1], 2002, 1), c(5, C[2], 2002, 1), c(5, C[3], 2002, 1),
  c(6, C[1], 2000, 1), c(6, C[2], 2001, 1), c(6, C[3], 2001, 1),
  c(6, C[4], 2002, 1), c(6, C[5], 2002, 2),
  c(6, C[1], 2001, 3)   # same species-cell as an earlier row, other year: counted once
))

test_that("taxon subsets include reptile classes and exclude missing class", {
  df <- data.frame(class = c("Amphibia", "Squamata", "Testudines", "Crocodylia",
                             "Aves", NA, "Insecta"),
                   kingdom = c(rep("Animalia", 6), "Animalia"))
  expect_equal(nrow(subset_taxon(df, "herps")), 4)
  expect_equal(nrow(subset_taxon(df, "aves")), 1)
  expect_equal(nrow(subset_taxon(df, "full")), 7)
  expect_false(anyNA(subset_taxon(df, "aves")$class))
  df2 <- data.frame(class = c("Magnoliopsida", NA), kingdom = c("Plantae", "Plantae"))
  expect_equal(nrow(subset_taxon(df2, "plants")), 2)
  expect_error(subset_taxon(df, "fish"))
})

test_that("year filter keeps 1990 onwards and drops missing years", {
  df <- data.frame(year = c(1989, 1990, 2020, NA))
  expect_equal(filter_years(df)$year, c(1990, 2020))
})

test_that("incidence stats count species-cell pairs once across years", {
  st <- incidence_stats(cube1)
  expect_equal(unlist(st), c(S = 6, T = 5, Q1 = 2, Q2 = 2, U = 14))
})

test_that("Chao2 matches the hand formula and vegan::specpool", {
  st <- incidence_stats(cube1)
  # S + (T-1)/T * Q1^2/(2 Q2) = 6 + 0.8 * 4/4 = 6.8
  expect_equal(chao2_estimate(st), 6.8)
  expect_equal(chao2_completeness(st), 6 / 6.8)
  pairs <- unique(cube1[, c("mgrscellcode", "specieskey")])
  m <- (unclass(table(pairs$mgrscellcode, pairs$specieskey)) > 0) * 1
  expect_equal(chao2_estimate(st), unname(vegan::specpool(m)$chao))
})

test_that("Chao2 with Q2 = 0 uses the Q1(Q1-1)/2 form, as vegan does", {
  # s1,s2,s3 in one cell each (4 cells); s4 in 3 cells; s5 in all 4
  cube <- make_cube(list(c(1, C[1], 2000, 1), c(2, C[2], 2000, 1), c(3, C[3], 2000, 1),
                         c(4, C[1], 2000, 1), c(4, C[2], 2000, 1), c(4, C[3], 2000, 1),
                         c(5, C[1], 2000, 1), c(5, C[2], 2000, 1), c(5, C[3], 2000, 1),
                         c(5, C[4], 2000, 1)))
  st <- incidence_stats(cube)
  expect_equal(c(st$S, st$T, st$Q1, st$Q2), c(5, 4, 3, 0))
  # 5 + 3/4 * 3*2/2 = 7.25
  expect_equal(chao2_estimate(st), 7.25)
  pairs <- unique(cube[, c("mgrscellcode", "specieskey")])
  m <- (unclass(table(pairs$mgrscellcode, pairs$specieskey)) > 0) * 1
  expect_equal(chao2_estimate(st), unname(vegan::specpool(m)$chao))
})

test_that("Chao2 is undefined with a single occupied cell (no fallback)", {
  cube <- make_cube(list(c(1, C[1], 2000, 1), c(2, C[1], 2000, 1)))
  st <- incidence_stats(cube)
  expect_true(is.na(chao2_estimate(st)))
  expect_true(is.na(chao2_completeness(st)))
  expect_true(is.na(pooled_coverage(st)))
})

test_that("pooled coverage matches the Chao & Jost formula and iNEXT", {
  st <- incidence_stats(cube1)
  # 1 - (2/14) * [4*2 / (4*2 + 2*2)] = 1 - (1/7)*(8/12) = 1 - 2/21
  expect_equal(pooled_coverage(st), 1 - 2 / 21)
  pairs <- unique(cube1[, c("specieskey", "mgrscellcode")])
  m <- (unclass(table(pairs$specieskey, pairs$mgrscellcode)) > 0) * 1L
  sc <- suppressWarnings(iNEXT::DataInfo(list(a = m), datatype = "incidence_raw"))$SC
  expect_equal(pooled_coverage(st), as.numeric(sc), tolerance = 1e-4)  # iNEXT rounds to 4 dp
})

test_that("annual coverage table records cells and matrix size per year", {
  tab <- annual_coverage_table(cube1)
  expect_equal(tab$year, 2000:2002)
  # 2000: s1,s2,s3,s6 in cells C1,C2 -> 4 species x 2 cells
  expect_equal(tab$n_cells[1], 2); expect_equal(tab$n_species[1], 4)
  expect_equal(tab$matrix_size[1], 8)
})

test_that("single-cell years get SC = 1 and only the min_cells rule drops them", {
  one_cell <- make_cube(lapply(1:6, function(s) c(s, C[1], 2000, 1)))   # 6 spp, 1 cell
  two_cell <- make_cube(list(c(1, C[1], 2001, 1), c(2, C[1], 2001, 1), c(3, C[1], 2001, 1),
                             c(4, C[2], 2001, 1), c(5, C[2], 2001, 1), c(6, C[2], 2001, 1)))
  tab <- annual_coverage_table(rbind(one_cell, two_cell))
  expect_equal(tab$sc[tab$year == 2000], 1)
  sc01 <- tab$sc[tab$year == 2001]
  expect_equal(mean_annual_coverage(tab, "b3gbi"), mean(c(1, sc01)))
  expect_equal(mean_annual_coverage(tab, "min_cells"), sc01)
})

test_that("mean annual coverage under the b3gbi rule equals b3gbi::completeness_ts", {
  skip_if_not_installed("b3gbi")
  # completeness_ts() reads lazy-loaded package data, so b3gbi must be attached.
  suppressPackageStartupMessages(library(b3gbi))
  set.seed(42)
  cells <- sprintf("32UNU%03d%03d", rep(1:6, 6), rep(1:6, each = 6))
  rows <- expand.grid(specieskey = 1:25, mgrscellcode = cells, year = 2000:2004,
                      stringsAsFactors = FALSE)
  rows <- rows[runif(nrow(rows)) < 0.12, ]
  rows$occurrences <- rpois(nrow(rows), 2) + 1
  rows$species <- paste0("sp", rows$specieskey)
  rows$kingdom <- "Animalia"; rows$class <- "Aves"
  rows$mincoordinateuncertaintyinmeters <- 50
  ours <- mean_annual_coverage(annual_coverage_table(rows), "b3gbi")
  cube <- suppressMessages(suppressWarnings(b3gbi::process_cube(rows,
            grid_type = "mgrs", cols_cellCode = "mgrscellcode", cols_speciesKey = "specieskey",
            cols_scientificName = "species", cols_occurrences = "occurrences")))
  theirs <- mean(as.data.frame(suppressMessages(suppressWarnings(
              b3gbi::completeness_ts(cube)))[["data"]])$diversity_val, na.rm = TRUE)
  expect_equal(ours, theirs, tolerance = 1e-9)
})

test_that("annual series: total occurrences and distinct species per year", {
  s <- annual_series(cube1)
  expect_equal(s$year, 2000:2002)
  # 2000 rows: s1 1, s2 1, s3 2, s6 1 -> 5 occ, 4 spp
  # 2001 rows: s3 1, s4 1, s6 1, s6 1, s6 3 -> 7 occ, 3 spp
  # 2002 rows: s4 4, s5 1,1,1, s6 1,2 -> 10 occ, 3 spp
  expect_equal(s$total_occ, c(5, 7, 10))
  expect_equal(s$richness, c(4, 3, 3))
})

test_that("Cha distance: identical shapes 0, reversed shapes worked example", {
  expect_equal(cha_distance(c(1, 2, 3), c(10, 20, 30)), 0)
  # P = (0, .5, 1), Q = (1, .5, 0): num = 1 + 0 + 1 = 2
  # den = 1.25 + 1.25 - 0.25 = 2.25  ->  2 / 2.25
  expect_equal(cha_distance(c(1, 2, 3), c(3, 2, 1)), 2 / 2.25)
  expect_true(is.na(cha_distance(c(1, 2), c(2, 1))))          # < 3 points
  expect_true(is.na(cha_distance(c(5, 5, 5), c(1, 2, 3))))    # constant series
  set.seed(1); d <- replicate(200, cha_distance(runif(8), runif(8)))
  expect_true(all(d >= 0 & d <= 1))
})

test_that("binary Jaccard reproduces the D6.1 worker's statistic", {
  # D6.1 calculate_jaccard, copied verbatim for comparison
  d61 <- function(occ_vec, rich_vec) {
    if (length(occ_vec) < 3) return(NA_real_)
    occ_norm <- (occ_vec - min(occ_vec)) / (max(occ_vec) - min(occ_vec) + 1e-10)
    rich_norm <- (rich_vec - min(rich_vec)) / (max(rich_vec) - min(rich_vec) + 1e-10)
    occ_set <- which(occ_norm > mean(occ_norm)); rich_set <- which(rich_norm > mean(rich_norm))
    u <- length(union(occ_set, rich_set)); if (u == 0) return(NA_real_)
    1 - length(intersect(occ_set, rich_set)) / u
  }
  set.seed(7)
  for (i in 1:100) { p <- rpois(12, 20); q <- rpois(12, 8)
    expect_equal(binary_jaccard_distance(p, q), d61(p, q)) }
})

test_that("elbow threshold matches the D6.1 function", {
  d61 <- function(distances) {
    dist_clean <- sort(distances[!is.na(distances)]); n <- length(dist_clean)
    x <- (1:n) / n; y <- (dist_clean - min(dist_clean)) / (max(dist_clean) - min(dist_clean))
    dist_clean[which.max(abs(x - y) / sqrt(2))]
  }
  set.seed(3); d <- c(rbeta(500, 2, 5), NA, NA)
  expect_equal(elbow_threshold(d), d61(d))
  expect_error(elbow_threshold(1:5))
})

test_that("precision ratio: occurrence-weighted mean of areas, missing share", {
  df <- data.frame(mincoordinateuncertaintyinmeters = c(100, 1000, NA),
                   occurrences = c(3, 1, 6))
  # areas: pi*0.1^2 = 0.0314159, pi*1^2 = 3.14159; weights 3 and 1
  expected <- (3 * pi * 0.01 + 1 * pi * 1) / 4
  r <- precision_ratio(df, area_km2 = 2)
  expect_equal(r$mean_area_km2, expected)
  expect_equal(r$ratio, expected / 2)
  expect_equal(r$frac_occ_missing, 0.6)
  expect_true(is.na(precision_ratio(df[3, ], 2)$ratio))
})

test_that("D6.1 precision definition: unweighted mean radius, then area", {
  df <- data.frame(mincoordinateuncertaintyinmeters = c(100, 1000, NA))
  expect_equal(precision_ratio_d61(df, 1), pi * 0.55^2)
})

test_that("density is occurrences per km2", {
  df <- data.frame(occurrences = c(3, 7, NA))
  expect_equal(density_km2(df, 40), 0.25)
  expect_true(is.na(density_km2(df, NA)))
  expect_true(is.na(density_km2(df, 0)))
})

test_that("WKT lookup does not confuse site_34 with site_342", {
  d <- file.path(tempdir(), "wkt_test"); unlink(d, recursive = TRUE)
  dir.create(file.path(d, "South Africa"), recursive = TRUE)
  file.create(file.path(d, "South Africa", c("site_34_De_Hoop_Vlei.wkt", "site_342_De_Mond.wkt")))
  expect_equal(basename(find_wkt(d, "site_34")), "site_34_De_Hoop_Vlei.wkt")
  expect_equal(basename(find_wkt(d, "site_342")), "site_342_De_Mond.wkt")
  expect_length(find_wkt(d, "site_3"), 0)
})

test_that("site area is geodesic km2", {
  f <- tempfile(fileext = ".wkt")
  writeLines("POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))", f)
  # Expected values from the closed-form WGS84 area of a 1 x 1 degree cell,
  # A = dlon * b^2 * [q(lat2) - q(lat1)] / 2, computed independently in Python.
  # That formula assumes edges along parallels; sf treats edges as geodesics,
  # which bow slightly poleward, so agreement is ~3e-5 rather than exact.
  expect_equal(site_area_km2(f), 12308.46, tolerance = 1e-4)
  writeLines("POLYGON ((10 60, 11 60, 11 61, 10 61, 10 60))", f)
  expect_equal(site_area_km2(f), 6123.14, tolerance = 1e-4)
})

test_that("species are identified by name, so two keys for one species count once", {
  cube <- make_cube(list(c(1, C[1], 2020, 1), c(1, C[2], 2020, 1), c(2, C[3], 2024, 1)))
  cube$species <- c("Palpopleura lucia", "Palpopleura lucia", "Palpopleura lucia")
  cube$specieskey <- c("1428955", "1428955", "4C7CK")   # backbone key, then new key
  expect_equal(incidence_stats(cube)$S, 1)
  expect_equal(incidence_stats(cube, "specieskey")$S, 2)
  expect_equal(names_with_multiple_keys(cube), 1)
  expect_equal(annual_series(cube)$richness, c(1, 1))
})

test_that("rows without a species name are excluded from incidence counts", {
  cube <- make_cube(list(c(1, C[1], 2020, 1), c(2, C[2], 2020, 1)))
  cube$species[2] <- ""
  expect_equal(incidence_stats(cube)$S, 1)
})

test_that("key harmonisation gives one key per name, preferring numeric keys", {
  df <- data.frame(species = c("A a", "A a", "A a", "B b", "B b", ""),
                   specieskey = c("4C7CK", "1428955", "4C7CK", "ZKQD", "ZKQD", "999"),
                   stringsAsFactors = FALSE)
  h <- harmonise_species_keys(df)
  expect_equal(h$specieskey, c("1428955", "1428955", "1428955", "ZKQD", "ZKQD", "999"))
  expect_equal(names_with_multiple_keys(h[h$species != "", ]), 0)
})

test_that("after harmonisation b3gbi counts a two-key species once", {
  skip_if_not_installed("b3gbi")
  suppressPackageStartupMessages(library(b3gbi))
  d <- data.frame(year = c(2023, 2023, 2024), species = "Palpopleura lucia",
                  specieskey = c("1428955", "4C7CK", "4C7CK"),
                  mgrscellcode = c("29PNM123906", "29PNM124907", "29PNM123906"),
                  occurrences = 1, stringsAsFactors = FALSE)
  rich <- function(x) {
    cube <- suppressMessages(suppressWarnings(b3gbi::process_cube(x, cols_cellCode = "mgrscellcode",
              cols_speciesKey = "specieskey", cols_scientificName = "species")))
    as.data.frame(b3gbi::obs_richness_ts(cube)$data)$diversity_val
  }
  expect_equal(rich(d), c(2, 1))                           # the problem
  expect_equal(rich(harmonise_species_keys(d)), c(1, 1))   # fixed
})
