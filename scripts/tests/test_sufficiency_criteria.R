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
  expect_error(subset_taxon(df, "fish"))            # no phylum column
  expect_error(subset_taxon(df, "trees"))
})

test_that("fish include chordates without a class (old backbone) but no other groups", {
  df <- data.frame(phylum = c("Chordata", "Chordata", "Chordata", "Chordata", "Chordata", "Arthropoda", "Chordata", "Chordata", "Chordata"),
                   class = c("", NA, "Elasmobranchii", "Actinopterygii", "Aves", "", "Ascidiacea", "Petromyzonti", ""),
                   order = c("Perciformes", "Salmoniformes", "Rajiformes", "Perciformes", "Anseriformes", "Araneae", "Stolidobranchia", "Petromyzontiformes", "Copelata"),
                   kingdom = "Animalia")
  expect_equal(which(rownames(df) %in% rownames(subset_taxon(df, "fish"))), c(1, 2, 3, 4, 8))
  expect_error(subset_taxon(df[, c("phylum", "class", "kingdom")], "fish"))   # needs order
  # the classes of the new backbone in the fish download
  nb <- data.frame(phylum = "Chordata", order = "x", kingdom = "Animalia",
                   class = c("Teleostei", "Elasmobranchii", "Chondrostei", "Holostei", "Holocephali", "Petromyzonti",
                             "Myxini", "Dipneusti", "Cladistii", "Aves", "Mammalia", "Reptilia", "Amphibia"))
  expect_equal(subset_taxon(nb, "fish")$class, nb$class[1:9])
})

test_that("full is no longer screened by default but can still be asked for", {
  expect_false("full" %in% SUBSETS)
  expect_setequal(SUBSETS, c("aves", "mammalia", "herps", "insects", "plants", "fish"))
  expect_true(all(c("full", SUBSETS) %in% ALL_SUBSETS))
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
  # incidences = species x cell pairs: 2000 has s1:C1, s2:C2, s3:C1, s6:C1 -> 4
  expect_equal(tab$n_incidences[1], 4)
})

test_that("the informative rule keeps only years with enough cells and incidences", {
  tab <- data.frame(year = 1:4, n_cells = c(19, 20, 40, 20), n_incidences = c(500, 49, 50, 50),
                    matrix_size = 1e3, sc = c(1, 0.99, 0.80, 0.90))
  expect_equal(informative_years(tab), c(FALSE, FALSE, TRUE, TRUE))
  expect_equal(mean_annual_coverage(tab, "informative"), mean(c(0.80, 0.90)))
  expect_true(is.na(mean_annual_coverage(tab[1:2, ], "informative")))   # no informative year: undefined
  expect_equal(informative_years(tab, min_cells = 0, min_incidences = 0), rep(TRUE, 4))
})

test_that("annual coverage equals iNEXT::DataInfo on the dense species x cell matrix", {
  dense_sc <- function(y) {
    pairs <- unique(y[, c("species", "mgrscellcode")])
    m <- (unclass(table(pairs$species, pairs$mgrscellcode)) > 0) * 1L
    tryCatch(as.numeric(suppressWarnings(iNEXT::DataInfo(list(a = m), datatype = "incidence_raw"))$SC),
             error = function(e) NA_real_)
  }
  set.seed(11)
  n_checked <- 0
  for (k in 1:300) {
    n_sp <- sample(1:30, 1); n_cell <- sample(1:12, 1); p <- runif(1, 0.02, 0.6)
    g <- expand.grid(species = paste0("sp", 1:n_sp), mgrscellcode = paste0("C", 1:n_cell),
                     stringsAsFactors = FALSE)
    g <- g[runif(nrow(g)) < p, ]
    if (nrow(g) == 0) next
    g$year <- 2000; g$occurrences <- 1
    tab <- annual_coverage_table(g)
    expect_equal(tab$sc, dense_sc(g))
    expect_equal(tab$matrix_size, length(unique(g$species)) * length(unique(g$mgrscellcode)))
    n_checked <- n_checked + 1
  }
  expect_gt(n_checked, 250)
  # edge cases: 1 species x 1 cell, 1 species x many cells, many species x 1 cell, all singletons
  e <- list(data.frame(species = "a", mgrscellcode = "C1"),
            data.frame(species = "a", mgrscellcode = paste0("C", 1:5)),
            data.frame(species = paste0("s", 1:5), mgrscellcode = "C1"),
            data.frame(species = paste0("s", 1:4), mgrscellcode = paste0("C", 1:4)),
            data.frame(species = c("a", "a", "b", "b"), mgrscellcode = c("C1", "C2", "C1", "C2")))
  for (g in e) {
    g$year <- 2000; g$occurrences <- 1
    expect_equal(annual_coverage_table(g)$sc, dense_sc(g))
  }
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
  # weighted median radius: weights 3 (100 m) and 1 (1000 m) -> 100 m
  expect_equal(r$median_area_km2, pi * 0.1^2)
  expect_equal(r$ratio_median, pi * 0.1^2 / 2)
  expect_true(is.na(precision_ratio(df[3, ], 2)$ratio_median))
})

test_that("weighted median", {
  expect_equal(weighted_median(c(3, 1, 2), c(1, 1, 1)), 2)
  expect_equal(weighted_median(c(1, 2, 3, 4), c(1, 1, 1, 1)), 2)   # lower median
  expect_equal(weighted_median(c(1, 10), c(1, 5)), 10)
  expect_equal(weighted_median(c(5, 1, 10), c(2, 1, 1)), 5)
  set.seed(9); x <- runif(101)
  expect_equal(weighted_median(x, rep(1, 101)), median(x))
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

test_that("harmonisation merges rows that collide, so b3gbi keeps every occurrence", {
  skip_if_not_installed("b3gbi")
  suppressPackageStartupMessages(library(b3gbi))
  # same cell and year under two keys; b3gbi drops identical rows after harmonising
  d <- data.frame(year = c(2019, 2019, 2019, 2020), species = c("P l", "P l", "Q q", "P l"),
                  specieskey = c("11345536", "1579202", "7", "1579202"),
                  mgrscellcode = c("30UWB024930", "30UWB024930", "30UWB024930", "30UWB024930"),
                  occurrences = c(1, 1, 3, 2), mincoordinateuncertaintyinmeters = c(70.7, 2, NA, 5),
                  stringsAsFactors = FALSE)
  h <- harmonise_species_keys(d)
  expect_equal(nrow(h), 3)
  expect_equal(h$occurrences, c(2, 3, 2))
  expect_equal(h$mincoordinateuncertaintyinmeters, c(2, NA, 5))
  expect_equal(sum(h$occurrences), sum(d$occurrences))
  occ <- function(x) {
    cube <- suppressMessages(suppressWarnings(b3gbi::process_cube(x, cols_cellCode = "mgrscellcode",
              cols_speciesKey = "specieskey", cols_scientificName = "species")))
    as.data.frame(b3gbi::total_occ_ts(cube)$data)$diversity_val
  }
  expect_equal(occ(h), c(5, 2))
  # without duplicates the data are returned unchanged
  expect_identical(merge_duplicate_cube_rows(d[3:4, ]), d[3:4, ])
})

test_that("MAGV reliability: rolling mean |change| and sd over 5 years", {
  magv <- c("scripts/lib/magv.R", "../lib/magv.R")
  source(magv[file.exists(magv)][1])
  n <- c(10, 20, 10, 10, 15, 15, 15, 15, 15, 15, 16)
  s <- magv_series(1990:2000, n)
  rc <- c(NA, 1, -0.5, 0, 0.5, 0, 0, 0, 0, 0, 1 / 15)
  expect_equal(s$relative_change, rc)
  expect_true(all(is.na(s$reliability[1:4])))
  w <- rc[1:5]; expect_equal(s$lma_abs[5], mean(abs(w), na.rm = TRUE))
  expect_equal(s$volatility[5], sd(w, na.rm = TRUE))
  w <- rc[6:10]; expect_equal(s$reliability[10], 1)          # flat window
  expect_equal(s$phase[10], "Plateau of Stability")
  expect_equal(s$phase[1], "Discovery Phase")
  # unsorted input is sorted by year
  expect_equal(magv_series(rev(1990:2000), rev(n))$reliability, s$reliability)
  m <- magv_summary(s)
  st <- which(s$phase == "Plateau of Stability")
  expect_equal(m$stable_run_length, max(rle(s$phase == "Plateau of Stability")$lengths[rle(s$phase == "Plateau of Stability")$values]))
  expect_equal(m$has_usable_window, m$stable_run_length >= 5)
  expect_equal(m$n_years, 11)
  flat <- magv_summary(magv_series(2000:2003, c(1, 2, 3, 4)))   # shorter than the window
  expect_equal(flat$stable_run_length, 0); expect_true(is.na(flat$mean_reliability))
})

test_that("site cube discovery: one file per site, the largest, first on ties", {
  root <- file.path(tempdir(), "disc"); unlink(root, recursive = TRUE)
  d1 <- file.path(root, "ramsar_site_data_100m_europe", "France"); dir.create(d1, recursive = TRUE)
  d2 <- file.path(root, "ramsar_site_data_100m_africa", "Chad"); dir.create(d2, recursive = TRUE)
  w <- function(path, n) writeLines(rep("x", n), path)
  w(file.path(d1, "site_5_A_data.csv"), 10); w(file.path(d1, "site_5_B_data.csv"), 50)
  w(file.path(d1, "site_50_C_data.csv"), 5); w(file.path(d1, "site_5_notes.txt"), 500)
  w(file.path(d2, "site_7_X_data.csv"), 20); w(file.path(d1, "site_7_Y_data.csv"), 20)
  x <- discover_site_cubes(root)
  expect_setequal(x$site_id, c("site_5", "site_50", "site_7"))
  expect_equal(basename(x$file_path[x$site_id == "site_5"]), "site_5_B_data.csv")
  expect_equal(basename(x$file_path[x$site_id == "site_7"]), "site_7_X_data.csv")  # africa comes first
  expect_equal(x$cont[x$site_id == "site_7"], "africa")
})

# ---- v1 indicators (scripts/lib/indicators_v1.R) ----------------------------
ind_lib <- c("scripts/lib/indicators_v1.R", "../lib/indicators_v1.R")
source(ind_lib[file.exists(ind_lib)][1])

test_that("MGRS 100 m codes map to the containing 1 km cell", {
  expect_equal(mgrs_to_1km(c("32VNM641645", "4QFJ123678", "32VNM649649")),
               c("32VNM6464", "4QFJ1267", "32VNM6464"))
  expect_error(mgrs_to_1km(c("32VNM6464")))
})

test_that("recode_grid leaves 100 m cubes unchanged and merges cells on the 1 km grid", {
  expect_identical(recode_grid(cube1, "100m"), cube1)
  k <- recode_grid(cube1, "1km")
  expect_equal(unique(k$mgrscellcode), "32UNU0000")   # all five 100 m cells lie in one 1 km cell
  expect_equal(k[, setdiff(names(k), "mgrscellcode")], cube1[, setdiff(names(cube1), "mgrscellcode")])
  st <- incidence_stats(k)
  expect_equal(c(st$S, st$T), c(6, 1))
  expect_true(is.na(chao2_completeness(st)))           # one sampling unit: undefined
  blank <- cube1[1:2, ]; blank$mgrscellcode[1] <- ""
  expect_equal(recode_grid(blank, "1km")$mgrscellcode, c("", "32UNU0000"))
})

test_that("simple indicators: occurrences, richness and cumulative richness by name", {
  d <- data.frame(year = c(2000, 2000, 2000, 2001, 2001, 2002),
                  species = c("A a", "B b", "A a", "A a", "C c", "B b"),
                  specieskey = c("1", "2", "1", "9", "3", "2"),
                  mgrscellcode = c("32VNM641645", "32VNM641645", "32VNM642646", "32VNM700700", "32VNM641645", "32VNM641645"),
                  occurrences = c(1, 2, 3, 4, 5, 6), stringsAsFactors = FALSE)
  s <- simple_indicators(d, "100m")
  expect_equal(s$total_occ, c(6, 9, 6))
  expect_equal(s$obs_richness, c(2, 2, 1))        # "A a" counted once despite two keys
  expect_equal(s$cum_richness, c(2, 3, 3))
  expect_equal(s$n_cells, c(2, 2, 1))
  s1 <- simple_indicators(d, "1km")
  expect_equal(s1$n_cells, c(1, 2, 1))             # 641645 and 642646 share 1 km cell 6464
  expect_equal(s1$total_occ, s$total_occ)
})

test_that("Hill numbers from frequencies equal b3gbi's hill*_ts on the 100 m grid", {
  skip_if_not_installed("b3gbi")
  suppressPackageStartupMessages(library(b3gbi))
  set.seed(7)
  cells <- sprintf("31TFJ%03d%03d", rep(1:8, 8) + 100, rep(1:8, each = 8) + 400)
  rows <- expand.grid(specieskey = 1:30, mgrscellcode = cells, year = 2000:2004, stringsAsFactors = FALSE)
  rows <- rows[runif(nrow(rows)) < 0.08, ]
  rows$occurrences <- rpois(nrow(rows), 2) + 1
  rows$species <- paste0("sp", rows$specieskey); rows$specieskey <- as.character(rows$specieskey)
  rows$kingdom <- "Animalia"; rows$class <- "Aves"; rows$mincoordinateuncertaintyinmeters <- 50
  cube <- suppressMessages(suppressWarnings(b3gbi::process_cube(rows, grid_type = "mgrs",
            cols_cellCode = "mgrscellcode", cols_speciesKey = "specieskey", cols_scientificName = "species",
            cols_occurrences = "occurrences")))
  ours <- hill_indicators(rows, "100m", level = 0.95)
  for (qq in 0:2) {
    b <- as.data.frame(suppressMessages(suppressWarnings(do.call(paste0("hill", qq, "_ts"),
           list(data = cube, ci_type = "none", cell_size = "grid"))))$data)
    b <- b[!is.na(b$diversity_val), ]
    o <- ours[ours$q == qq, ]
    expect_equal(o$year, b$year)
    expect_equal(o$qD, b$diversity_val, tolerance = 1e-10)
  }
  # intervals: point estimates unchanged, bounds present where the bootstrap succeeds
  set.seed(1); ci <- hill_indicators(rows, "100m", level = 0.95, nboot = 20)
  expect_equal(ci$qD, ours$qD)
  expect_true(all(ci$qD.LCL[ci$ci_status == "ok"] <= ci$qD[ci$ci_status == "ok"] + 1e-9))
  expect_true(all(is.na(ci$qD.LCL[ci$ci_status == "failed"])))
  # selected years: identical to the same years of a run over all years, intervals included
  part <- hill_indicators(rows, "100m", level = 0.95, nboot = 20, years = 2003:2004)
  expect_equal(sort(unique(part$year)), c(2003, 2004))
  expect_equal(part, ci[ci$year %in% 2003:2004, ], ignore_attr = TRUE)
  expect_equal(nrow(hill_indicators(rows, "100m", level = 0.95, years = 1990)), 0)
})
