# Unit tests for scripts/lib/cube_update.R and scripts/lib/taxonomy_match.R
# Run from the repository root:
#   Rscript -e 'testthat::test_file("scripts/tests/test_cube_update.R")'
#
# The API responses in tests/fixtures/ are real responses of the GBIF v2
# species-match API (new backbone), shortened to the fields used.

library(testthat)
suppressPackageStartupMessages({ library(data.table); library(sf); library(mgrs) })
here <- if (file.exists("scripts/lib/cube_update.R")) "scripts" else ".."
source(file.path(here, "lib/cube_update.R"))
source(file.path(here, "lib/taxonomy_match.R"))
source(file.path(here, "lib/sufficiency_criteria.R"))
fixture <- function(f) paste(readLines(file.path(here, "tests/fixtures", f), warn = FALSE), collapse = "")

cube_row <- function(species, key, year, cell, occ = 1, unc = 1000, kingdom = "Animalia", cls = "Aves") {
  x <- as.list(setNames(rep("1", length(CUBE_COLS)), CUBE_COLS))
  x$kingdom <- kingdom; x$class <- cls; x$species <- species; x$specieskey <- key; x$year <- as.character(year)
  x$mgrscellcode <- cell; x$occurrences <- as.character(occ); x$mincoordinateuncertaintyinmeters <- as.character(unc)
  as.data.table(x)
}

test_that("site ids are read without splitting on dots", {
  expect_equal(site_id_of("site_1979_Autertal_-_St._Lorenzener_Hochmoor.wkt"), "site_1979")
  expect_equal(site_id_of(c("a/site_828_Lagoa_de_St._Andre_data.csv", "x/other.csv")), c("site_828", NA))
})

test_that("outline lookup does not confuse site_1 with site_10 and prefers the cube's country", {
  w <- c("wkt/A/site_10_X.wkt", "wkt/A/site_1_Y.wkt", "wkt/B/site_1_Y.wkt")
  expect_equal(find_site_wkt("site_1", w), "wkt/A/site_1_Y.wkt")
  expect_equal(find_site_wkt("site_1", w, country = "B"), "wkt/B/site_1_Y.wkt")
  expect_equal(find_site_wkt("site_10", w), "wkt/A/site_10_X.wkt")
  expect_true(is.na(find_site_wkt("site_2", w)))
})

test_that("needs_patch tells old-backbone and updated cubes apart", {
  f <- tempfile(fileext = ".csv")
  fwrite(rbind(cube_row("A b", "2480260", 2023, "30UVB883951"), cube_row("A b", "2480260", 2024, "30UVB883951")), f)
  expect_true(needs_patch(f))                         # 2024 rows with numeric keys
  fwrite(rbind(cube_row("A b", "2480260", 2023, "30UVB883951"), cube_row("A b", "7CYVX", 2024, "30UVB883951"),
               cube_row("C d", "69737", 2024, "30UVB883951"), cube_row("E f", "3J8XQ", 2025, "30UVB883951")), f)
  expect_false(needs_patch(f))                        # mostly new keys; one happens to be all digits
  fwrite(cube_row("A b", "2480260", 2020, "30UVB883951"), f)
  expect_true(needs_patch(f))                         # nothing from 2024 on
})

test_that("cells are assigned to the outline that contains them", {
  inside <- latlng_to_mgrs(51.40, -3.10, precision = 3)
  outside <- latlng_to_mgrs(51.80, -3.10, precision = 3)
  sq <- st_sf(site_id = c("site_1", "site_2"),
              geometry = st_sfc(st_polygon(list(rbind(c(-3.2, 51.3), c(-3.0, 51.3), c(-3.0, 51.5), c(-3.2, 51.5), c(-3.2, 51.3)))),
                                st_polygon(list(rbind(c(-3.15, 51.35), c(-3.05, 51.35), c(-3.05, 51.45), c(-3.15, 51.45), c(-3.15, 51.35)))),
                                crs = 4326))
  a <- assign_cells(c(inside, outside, "not_a_cell"), sq)
  expect_setequal(a$site_id, c("site_1", "site_2"))   # overlapping outlines: the cell belongs to both
  expect_true(all(a$mgrscellcode == inside))
  expect_true(substr(inside, 1, 5) %in% site_squares(st_geometry(sq)[1]))
})

test_that("stream_update reads past invalid bytes, filters by square and drops broken lines", {
  f <- tempfile(fileext = ".tsv")
  ln <- function(cell, sp) paste(c("Animalia", "N", "Chordata", "CH2", "Aves", "V2", "O", "o", "F", "f", "G", "g",
                                   sp, "K1", "2024", cell, "1", "1", "1", "1", "1", "1", "2", "60", "100"), collapse = "\t")
  con <- file(f, "wb")
  writeLines(paste(CUBE_COLS, collapse = "\t"), con, useBytes = TRUE)
  writeLines(ln("30UVB883951", "Keep one"), con, useBytes = TRUE)
  writeBin(c(charToRaw(ln("30UVB883952", "Bad ")), as.raw(c(0xff, 0xfe)), charToRaw("\n")), con)   # invalid UTF-8
  writeLines(ln("31UDQ000000", "Other square"), con, useBytes = TRUE)
  writeLines(ln("30UVB883953", "Keep two"), con, useBytes = TRUE)
  writeLines(paste(c("Animalia", "N", "Chordata", "CH2", "Aves", "V2", "O", "o", "F", "f", "G", "g", "Cut", "K", "2024", "30UVB883954", "1"), collapse = "\t"), con, useBytes = TRUE)  # truncated
  close(con)
  got <- list()
  n <- stream_update(f, "30UVB", on_chunk = function(dt) got[[length(got) + 1]] <<- dt, chunk_lines = 2, progress = FALSE)
  got <- rbindlist(got)
  expect_equal(n, 5)
  expect_setequal(got$mgrscellcode, c("30UVB883951", "30UVB883952", "30UVB883953"))
  expect_true("Keep two" %in% got$species)             # read on after the invalid bytes
})

test_that("patched_rows keeps earlier years and replaces the rest", {
  old <- rbind(cube_row("A b", "1", 2023, "30UVB883951"), cube_row("A b", "1", 2024, "30UVB883951"))
  old$ramsar_site_id <- "site_67"                      # extra site columns are dropped
  new <- rbind(cube_row("A b", "X", 2024, "30UVB883951", occ = 5), cube_row("C d", "Y", 2025, "30UVB883952"))
  r <- patched_rows(old, new, 2024)
  expect_equal(names(r), CUBE_COLS)
  expect_equal(r$year, c("2023", "2024", "2025"))
  expect_equal(r$occurrences[2], "5")
  expect_equal(row_keys(r), row_keys(new))
})

test_that("parse_match resolves synonyms and rejects higher-rank matches", {
  s <- parse_match(fixture("match_sylvia_curruca.json"))
  expect_equal(s$species, "Curruca curruca"); expect_equal(s$specieskey, "DDBMF")
  expect_equal(s$genus, "Curruca"); expect_equal(s$classkey, "V2"); expect_equal(s$status, "synonym")
  h <- parse_match(fixture("match_haematopus_ostralegus.json"))
  expect_equal(h$status, "accepted"); expect_equal(h$specieskey, "3J8XQ"); expect_equal(h$familykey, "625G4")
  p <- parse_match(fixture("match_ammophila_arenaria.json"))
  expect_equal(p$species, "Calamagrostis arenaria"); expect_equal(p$kingdom, "Plantae")
  x <- parse_match(fixture("match_xyzus_fakeus.json"))
  expect_equal(x$status, "unmatched"); expect_true(is.na(x$species))
  expect_equal(parse_match("not json")$status, "unmatched")
})

test_that("build_name_table uses the update data first, then the API, and resumes from its cache", {
  old <- data.table(kingdom = c("Animalia", "Animalia", "Animalia", "Plantae"), phylum = NA, class = NA,
                    species = c("Haematopus ostralegus", "Sylvia curruca", "Xyzus fakeus", "Ammophila arenaria"))
  new_rows <- as.data.table(parse_match(fixture("match_haematopus_ostralegus.json"))[, TAXON_COLS])
  calls <- 0
  fake <- function(urls) {
    calls <<- calls + length(urls)
    vapply(urls, function(u) {
      if (grepl("Sylvia", u)) fixture("match_sylvia_curruca.json")
      else if (grepl("Ammophila", u)) fixture("match_ammophila_arenaria.json")
      else if (grepl("Xyzus", u)) fixture("match_xyzus_fakeus.json") else NA_character_
    }, "", USE.NAMES = FALSE)
  }
  cache <- tempfile(fileext = ".csv")
  nt <- build_name_table(old, new_rows, cache, fetch = fake, progress = FALSE)
  expect_equal(calls, 3)                                # Haematopus came from the update data
  expect_equal(nt[old_species == "Sylvia curruca"]$species, "Curruca curruca")
  expect_equal(nt[old_species == "Haematopus ostralegus"]$status, "in update data")
  expect_equal(nt[old_species == "Xyzus fakeus"]$status, "unmatched")
  nt2 <- build_name_table(old, new_rows, cache, fetch = fake, progress = FALSE)
  expect_equal(calls, 3)                                # nothing asked again
  expect_equal(nrow(nt2), 4)
  expect_match(match_url("Sylvia curruca", "Animalia"), "scientificName=Sylvia%20curruca&kingdom=Animalia")
})

test_that("harmonise_rows renames old-backbone rows and merges the duplicates this creates", {
  nt <- data.table(old_kingdom = "Animalia", old_species = c("Sylvia curruca", "Sylvia althaea", "Xyzus fakeus"),
                   rbind(parse_match(fixture("match_sylvia_curruca.json"))[, TAXON_COLS],
                         parse_match(fixture("match_sylvia_curruca.json"))[, TAXON_COLS],
                         parse_match(fixture("match_xyzus_fakeus.json"))[, TAXON_COLS]),
                   status = c("synonym", "synonym", "unmatched"))
  x <- rbind(cube_row("Sylvia curruca", "2493052", 2010, "30UVB883951", occ = 2, unc = 1000),
             cube_row("Sylvia althaea", "5231198", 2010, "30UVB883951", occ = 3, unc = 100),   # lumped into the same species
             cube_row("Sylvia curruca", "2493052", 2011, "30UVB883951", occ = 1),
             cube_row("Xyzus fakeus", "999", 2010, "30UVB883951"),
             cube_row("Curruca curruca", "DDBMF", 2024, "30UVB883951", occ = 4))
  h <- harmonise_rows(x, nt)
  r <- h$rows
  expect_equal(h$n_merged, 1)
  expect_equal(h$n_renamed, 3)
  y2010 <- r[year == "2010" & species == "Curruca curruca"]
  expect_equal(nrow(y2010), 1); expect_equal(y2010$occurrences, "5"); expect_equal(y2010$mincoordinateuncertaintyinmeters, "100")
  expect_equal(y2010$genus, "Curruca"); expect_equal(y2010$specieskey, "DDBMF")
  expect_equal(r[species == "Xyzus fakeus"]$specieskey, "999")   # unmatched names keep their old name and key
  expect_equal(sort(unique(r$species)), c("Curruca curruca", "Xyzus fakeus"))
  expect_equal(nrow(r), 4)
})

test_that("only names in the screened groups are matched", {
  old <- data.table(kingdom = c("Animalia", "Animalia", "Animalia", "Plantae", "Fungi", "Animalia", "Animalia"),
                    phylum = c("Chordata", "Chordata", "Chordata", "Tracheophyta", "Basidiomycota", "Mollusca", "Chordata"),
                    class = c("Aves", "", "Squamata", "Magnoliopsida", "Agaricomycetes", "Gastropoda", "Elasmobranchii"),
                    order = c("Anseriformes", "Salmoniformes", "Squamata", "Poales", "Agaricales", "Stylommatophora", "Rajiformes"),
                    species = c("A a", "Salmo trutta", "Natrix natrix", "P p", "F f", "M m", "Raja clavata"))
  expect_equal(old[in_subsets(old)]$species, c("A a", "Salmo trutta", "Natrix natrix", "P p", "Raja clavata"))
  expect_equal(old[in_subsets(old, c("aves", "plants"))]$species, c("A a", "P p"))
  expect_equal(nrow(old[in_subsets(old[0])]), 0)
})
