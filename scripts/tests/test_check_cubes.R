# Tests for scripts/check_site_cubes.R and scripts/lib/cube_check.R.
# Run from the repository root:
#   Rscript -e 'testthat::test_file("scripts/tests/test_check_cubes.R")'
#
# A small made-up world (three sites, two continental files, two update files)
# is built in a temporary folder, cubes with known errors are written, and the
# check must find exactly those errors.

library(testthat)
suppressPackageStartupMessages({ library(data.table); library(sf); library(mgrs); library(stringi) })
here <- if (file.exists("scripts/lib/cube_check.R")) "scripts" else ".."
source(file.path(here, "lib/sufficiency_criteria.R"))
source(file.path(here, "lib/cube_check.R"))

test_that("byte ranges read every line exactly once, whatever the range size", {
  f <- tempfile(fileext = ".tsv")
  body <- c(paste(CHK_COLS, collapse = "\t"),
            sprintf("A\t1\tC\t2\tAves\t3\to\t4\tf\t5\tg\t6\tSpécies %d\tK%d\t2020\t31UES%06d\t1\t1\t1\t1\t1\t1\t%d\t60\t100", 1:500, 1:500, 1:500, 1:500))
  writeLines(body, f, useBytes = TRUE)
  for (pb in c(37, 1000, 4096, file.size(f))) {
    plan <- chk_plan_parts(f, pb)
    got <- list(); parts <- list()
    for (i in seq_len(nrow(plan))) {
      r <- chk_read_range(f, plan$start[i], plan$end[i], on_chunk = function(l) got[[length(got) + 1]] <<- l, chunk_lines = 7)
      parts[[i]] <- data.table(part = i, first = r$first, after = r$after, size = plan$size[i], eol_bytes = 1L, final_eol = TRUE)
    }
    expect_identical(charToRaw(paste(unlist(got), collapse = "\n")), charToRaw(paste(body[-1], collapse = "\n")))
    expect_true(chk_check_tiling(rbindlist(parts))$ok)
  }
  bad <- rbindlist(parts); bad$after[1] <- bad$after[1] - 1
  expect_false(chk_check_tiling(bad)$ok)
})

test_that("lines with the wrong number of fields are counted, invalid UTF-8 is cleaned", {
  ok <- paste(c("Animalia", "1", "Chordata", "44", "Aves", "212", "O", "1", "F", "1", "G", "1", "Ab c", "9", "2020", "31UES907782", rep("1", 6), "2", "60", "15.0"), collapse = "\t")
  bad_utf <- paste0(sub("Ab c", "Ab ", ok, fixed = TRUE), "")
  bad_utf <- rawToChar(c(charToRaw(sub("Ab c", "Ab", ok, fixed = TRUE))))
  bad_utf <- sub("Ab", paste0("Ab", rawToChar(as.raw(0xff))), bad_utf, fixed = TRUE)
  p <- chk_parse_lines(c(ok, "", "too\tfew", bad_utf))
  expect_equal(nrow(p$rows), 2); expect_equal(p$n_empty, 1); expect_equal(p$n_wrong_fields, 1); expect_equal(p$n_bad_utf8, 1)
  expect_equal(p$rows$species, c("Ab c", "Ab"))
  expect_equal(chk_normalise(p$rows)$mincoordinateuncertaintyinmeters, c("15", "15"))
})

# ---- a made-up world ----------------------------------------------------------
world <- function() {
  d <- tempfile("cubecheck_"); dir.create(d)
  ex <- file.path(d, "inst/extdata"); dir.create(ex, recursive = TRUE)
  sq <- function(x0, y0, w) st_polygon(list(rbind(c(x0, y0), c(x0 + w, y0), c(x0 + w, y0 + w), c(x0, y0 + w), c(x0, y0))))
  # site 1: one part; site 2: two parts (the shapefile has one row per part); site 3: no cube
  shp <- st_sf(v_idris = 1:4, ramsarid = c(1L, 2L, 2L, 3L), officialna = c("One", "Two", "Two", "Three"),
               iso3 = "BEL", country_en = "Belgium", area_off = c(100, 50, 50, 10),
               geometry = st_sfc(sq(4.00, 51.00, 0.05), sq(4.20, 51.00, 0.02), sq(4.22, 51.00, 0.02), sq(4.40, 51.00, 0.02), crs = 4326))
  dir.create(file.path(ex, "Ramsar_boundaries_160626"))
  suppressWarnings(st_write(shp, file.path(ex, "Ramsar_boundaries_160626/features_publishedPolygon.shp"), quiet = TRUE))
  w <- file.path(ex, "ramsar_sites_wkt/Belgium"); dir.create(w, recursive = TRUE)
  writeLines(st_as_text(st_geometry(shp)[[1]]), file.path(w, "site_1_One.wkt"))
  writeLines(st_as_text(st_union(st_geometry(shp)[2:3])), file.path(w, "site_2_Two.wkt"))
  writeLines(st_as_text(sq(4.40, 51.00, 0.03)), file.path(w, "site_3_Three.wkt"))   # 0.03 instead of 0.02: should be flagged
  cell <- function(lng, lat) latlng_to_mgrs(lat, lng, precision = 3)
  row <- function(sp, key, cls, yr, c, occ = 1, kingdom = "Animalia", phylum = "Chordata")
    paste(c(kingdom, "1", phylum, "44", cls, "212", "O", "1", "F", "1", "G", "1", sp, key, yr, c, rep("1", 6), occ, "60", "15.0"), collapse = "\t")
  c1 <- cell(4.010, 51.010); c1b <- cell(4.030, 51.030); c2 <- cell(4.205, 51.005); c2b <- cell(4.235, 51.005)
  c3 <- cell(4.405, 51.005); out <- cell(4.150, 51.010)
  hdr <- paste(CHK_COLS, collapse = "\t")
  eu <- c(hdr,
          row("Anas a", "101", "Aves", 2014, c1), row("Anas a", "101", "Aves", 2015, c1, 3), row("Anas a", "101", "Aves", 2015, c1b),
          row("Natrix n", "202", "Squamata", 2016, c1b), row("Bufo b", "303", "Amphibia", 2016, c1),
          row("Anas a", "101", "Aves", 2016, c2), row("Anas a", "101", "Aves", 2017, c2b),
          row("Anas a", "101", "Aves", 2018, c3), row("Anas a", "101", "Aves", 2018, out),
          row("Anas a", "101", "Aves", 2024, c1),                                   # continental 2024: left out by design
          "broken\tline", row("Anas a", "101", "Aves", 2019, "31UZZ000000"),         # unreadable line; cell without coordinates
          row("Anas a", "101", "Aves", 2019, "notacell"))                           # invalid cell code
  as <- c(hdr, row("Bufo b", "303", "Amphibia", 2020, c1))                          # site 1 data filed under another continent
  up <- c(hdr, row("Anas a", "7CYVX", "Aves", 2024, c1, 2), row("Natrix n", "RT9", "Reptilia", 2025, c1b),
          row("Anas a", "7CYVX", "Aves", 2023, c1),                                 # update before 2024: left out
          row("Anas a", "7CYVX", "Aves", 2026, c1))                                 # after --to-year 2025: left out
  fish <- c(hdr, row("Salmo t", "FSH1", "Teleostei", 2024, c2),
            row("Anas a", "7CYVX", "Aves", 2024, c2))                               # not a fish: left out by --update-groups
  cd <- file.path(ex, "continental_gbif_data"); dir.create(cd)
  writeLines(eu, file.path(cd, "europe_100m_mgrs_all.csv")); writeLines(as, file.path(cd, "asia_100m_mgrs_all.csv"))
  writeLines(eu[1:3], file.path(cd, "asia_1km_mgrs_all.csv"))                         # not a 100 m cube: ignored
  ud <- file.path(ex, "ramsar_site_data_100m_world_2024_2026"); dir.create(ud)
  writeLines(up, file.path(ud, "update.csv")); writeLines(fish, file.path(ud, "fish.csv"))
  # expected cubes
  parse <- function(l) chk_parse_lines(l)$rows
  site1 <- parse(c(eu[2:6], as[2], up[2:3])); site2 <- parse(c(eu[7:8], fish[2]))
  cdir <- file.path(ex, "ramsar_site_data_100m_europe/Belgium"); dir.create(cdir, recursive = TRUE)
  # site 1 with errors: 2015 rows dropped, the snake dropped, one row outside the site added, one row doubled
  bad1 <- site1[!(year == "2015" | class %in% c("Squamata", "Reptilia"))]
  bad1 <- rbind(bad1, parse(eu[10]), bad1[1])
  fwrite(bad1, file.path(cdir, "site_1_One_data.csv"))
  # site 2 correct, written in the older format (site columns in front, numbers without ".0")
  s2 <- cbind(ramsar_site_id = "site_2_Two", site_name = "site_2_Two", site2)
  s2$mincoordinateuncertaintyinmeters <- "15"
  fwrite(s2, file.path(cdir, "site_2_Two_data.csv"))
  # a cube for a site that is not in the shapefile
  fwrite(site2[1], file.path(cdir, "site_99_Nowhere_data.csv"))
  list(dir = d, site1 = site1, site2 = site2, cdir = cdir)
}

run_check <- function(w, ...) {
  repo <- normalizePath(file.path(here, ".."))
  args <- c(file.path(repo, "scripts/check_site_cubes.R"), "--root", file.path(w$dir, "inst/extdata"), "--out", file.path(w$dir, "check"),
            "--update", paste(file.path(w$dir, "inst/extdata/ramsar_site_data_100m_world_2024_2026", c("update.csv", "fish.csv")), collapse = ","),
            "--update-groups", "all,fish", "--to-year", "2025",
            "--cores", "2", "--compare-cores", "1", "--part-gb", "3e-7", "--buckets", "3", ...)
  out <- withr::with_dir(repo, system2("Rscript", args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) stop(paste(out, collapse = "\n"))
  invisible(out)
}

test_that("the check finds every planted error and nothing else", {
  skip_if_not_installed("withr")
  w <- world()
  run_check(w, "--step", "all", "--tag", "t1")
  ck <- file.path(w$dir, "check")
  inv <- fread(file.path(ck, "inventory.csv"))
  expect_match(inv[site_id == "site_3"]$issues, "no cube file")
  expect_match(inv[site_id == "site_99"]$issues, "not in the shapefile")
  expect_equal(inv[site_id %in% c("site_1", "site_2")]$issues, c("", ""))
  oc <- fread(file.path(ck, "outline_check.csv"))
  expect_lt(oc[site_id == "site_1"]$rel_symdiff, 1e-9); expect_lt(oc[site_id == "site_2"]$rel_symdiff, 1e-9)
  expect_gt(oc[site_id == "site_3"]$rel_symdiff, 0.5)
  expect_gt(oc[site_id == "site_3"]$mean_shift_m, 400); expect_lt(oc[site_id == "site_3"]$mean_shift_m, 700)
  expect_lt(oc[site_id == "site_1"]$mean_shift_m, 1e-3)
  src <- fread(file.path(ck, "reference/sources.csv"))
  expect_true(all(src$read_completely))
  expect_equal(nrow(src), 4)                                         # asia_1km is not read
  expect_equal(src[basename(source) == "europe_100m_mgrs_all.csv"]$n_wrong_fields, 1)
  expect_equal(src[basename(source) == "europe_100m_mgrs_all.csv"]$n_invalid_code, 1)
  expect_equal(src[basename(source) == "europe_100m_mgrs_all.csv"]$n_no_coords, 1)
  expect_equal(src[basename(source) == "europe_100m_mgrs_all.csv"]$rows_other_period, 1)
  expect_equal(src[basename(source) == "update.csv"]$rows_other_period, 2)
  expect_equal(src[basename(source) == "fish.csv"]$rows_other_taxa, 1)
  expect_equal(src[basename(source) == "fish.csv"]$other_taxa_classes, "Aves:1")
  st <- fread(file.path(ck, "compare_t1/sites.csv"))
  s1 <- st[site_id == "site_1"]; s2 <- st[site_id == "site_2"]; s3 <- st[site_id == "site_3"]
  expect_equal(s1$rows_ref, nrow(w$site1))
  expect_equal(s1$missing_before, 3)                                 # two 2015 rows and the 2016 snake
  expect_equal(s1$missing_from, 1)                                   # the 2025 snake
  expect_equal(s1$extra_outside_site, 1)
  expect_equal(s1$extra_inside_site, 1)                              # the doubled row
  expect_equal(s1$dup_cell_year_key, 1)
  expect_equal(c(s2$rows_cube, s2$rows_ref, s2$rows_same), rep(nrow(w$site2), 3))
  expect_equal(s2$missing_before + s2$missing_from + s2$extra_inside_site + s2$extra_outside_site, 0)
  expect_equal(c(s3$n_cube_files, s3$rows_ref), c(0, 1))
  expect_false(st[site_id == "site_99"]$in_shapefile)
  yg <- fread(file.path(ck, "compare_t1/year_group.csv.gz"))
  expect_equal(yg[site_id == "site_1" & year == 2015 & group == "aves"]$rows_cube, 0)
  expect_equal(yg[site_id == "site_1" & year == 2015 & group == "aves"]$rows_ref, 2)
  expect_equal(yg[site_id == "site_2" & group == "fish"]$rows_ref, 1)
  rep <- readLines(file.path(ck, "REPORT_t1.md"))
  expect_true(any(grepl("Not ready for analysis", rep)))
  expect_true(any(grepl("cubes that differ from their sources (1)", rep, fixed = TRUE)))
  expect_true(any(grepl("Classes left out of fish.csv", rep, fixed = TRUE)))
  # a second run reuses everything and gives the same answer
  run_check(w, "--step", "compare", "--tag", "t1")
  setindex(st, NULL)
  expect_equal(fread(file.path(ck, "compare_t1/sites.csv")), st)
  # a cube that changes afterwards is compared again under the same tag
  Sys.sleep(1)
  fwrite(w$site2[-1], file.path(w$cdir, "site_2_Two_data.csv"))
  run_check(w, "--step", "compare", "--tag", "t1")
  s2 <- fread(file.path(ck, "compare_t1/sites.csv"))[site_id == "site_2"]
  expect_equal(c(s2$rows_cube, s2$missing_before), c(2, 1))
})

test_that("harmonised cubes are compared after renaming the reference with the name table", {
  skip_if_not_installed("withr")
  w <- world()
  run_check(w, "--step", "inventory"); run_check(w, "--step", "reference")
  # harmonise site 2 correctly: Anas a (old key 101) -> Mareca a (NEW1); both 2016 and 2017 rows keep their own cells
  nt <- data.table(old_kingdom = "Animalia", old_species = "Anas a", kingdom = "Animalia", kingdomkey = "N", phylum = "Chordata", phylumkey = "CH2",
                   class = "Aves", classkey = "V2", order = "O", orderkey = "1", family = "F", familykey = "1", genus = "Mareca",
                   genuskey = "1", species = "Mareca a", specieskey = "NEW1", status = "synonym")
  ntf <- file.path(w$dir, "name_table.csv"); fwrite(nt, ntf)
  h2 <- copy(w$site2); h2[specieskey == "101", `:=`(species = "Mareca a", specieskey = "NEW1", kingdomkey = "N", phylumkey = "CH2", classkey = "V2", genus = "Mareca")]
  fwrite(h2, file.path(w$cdir, "site_2_Two_data.csv"))
  run_check(w, "--step", "compare", "--tag", "h", "--name-table", ntf, "--sites", "site_2")
  s2 <- fread(file.path(w$dir, "check/compare_h/sites.csv"))
  expect_equal(c(s2$rows_cube, s2$rows_same), c(3, 3))
  # a wrong rename is found
  h2[specieskey == "NEW1", species := "Wrong name"]
  fwrite(h2, file.path(w$cdir, "site_2_Two_data.csv"))
  run_check(w, "--step", "compare", "--tag", "h2", "--name-table", ntf, "--sites", "site_2")
  s2 <- fread(file.path(w$dir, "check/compare_h2/sites.csv"))
  expect_equal(s2$extra_inside_site, 2); expect_equal(s2$missing_before, 2)
})

test_that("merging by the name table sums occurrences and keeps the smallest uncertainty", {
  ref <- chk_normalise(data.table(kingdom = "Animalia", kingdomkey = "1", phylum = "Chordata", phylumkey = "44", class = "Aves",
                                  classkey = "212", order = "O", orderkey = "1", family = "F", familykey = "1", genus = "Sylvia",
                                  genuskey = "1", species = c("Sylvia curruca", "Sylvia althaea"), specieskey = c("2493052", "5231198"),
                                  year = "2010", mgrscellcode = "30UVB883951", occurrences = c("2", "3"),
                                  mintemporaluncertainty = c("60", "86400"), mincoordinateuncertaintyinmeters = c("1000", "100.0")),
                       CHK_HARMONISED_COLS)
  nt <- data.table(old_kingdom = "Animalia", old_species = c("Sylvia curruca", "Sylvia althaea"), kingdom = "Animalia", kingdomkey = "N",
                   phylum = "Chordata", phylumkey = "CH2", class = "Aves", classkey = "V2", order = "O", orderkey = "PQ", family = "F",
                   familykey = "X", genus = "Curruca", genuskey = "G", species = "Curruca curruca", specieskey = "DDBMF", status = "synonym")
  h <- chk_apply_name_table(ref, nt)
  expect_equal(nrow(h), 1)
  expect_equal(h$occurrences, "5"); expect_equal(h$mincoordinateuncertaintyinmeters, "100"); expect_equal(h$mintemporaluncertainty, "60")
  expect_equal(h$species, "Curruca curruca")
})

test_that("cubes rebuilt by the patch script pass the check", {
  skip_if_not_installed("withr")
  w <- world()
  repo <- normalizePath(file.path(here, ".."))
  # cubes as the continental cubes made them (years before 2024), with one stale old-backbone 2024 row and a
  # 2026 row that the patch must remove
  s1 <- w$site1[as.numeric(year) < 2024]; s2 <- w$site2[as.numeric(year) < 2024]
  stale <- copy(s1[1]); stale[, year := "2024"]
  late <- copy(s1[1]); late[, year := "2026"]
  fwrite(rbind(s1, stale, late), file.path(w$cdir, "site_1_One_data.csv")); fwrite(s2, file.path(w$cdir, "site_2_Two_data.csv"))
  file.symlink(file.path(repo, "scripts"), file.path(w$dir, "scripts"))
  # a second copy of site 1's cube in the unknown folder (as for site 1979): moved away on install
  unk <- file.path(w$dir, "inst/extdata/ramsar_site_data_100m_unknown/x"); dir.create(unk, recursive = TRUE)
  file.copy(file.path(w$cdir, "site_1_One_data.csv"), file.path(unk, "site_1_One_data.csv"))
  orig <- tools::md5sum(file.path(w$cdir, c("site_1_One_data.csv", "site_2_Two_data.csv")))
  ud <- "inst/extdata/ramsar_site_data_100m_world_2024_2026"
  patch <- function(...) withr::with_dir(w$dir, system2("Rscript", c("scripts/patch_site_cubes_update.R", ...), stdout = TRUE, stderr = TRUE))
  out <- patch("--all", "--out-dir", "_local_data/cube_patch_all",
               "--update", paste(file.path(ud, c("update.csv", "fish.csv")), collapse = ","), "--update-groups", "all,fish",
               "--to-year", "2025", "--chunk-lines", "2")
  expect_null(attr(out, "status"))
  expect_equal(tools::md5sum(file.path(w$cdir, c("site_1_One_data.csv", "site_2_Two_data.csv"))), orig)   # cubes untouched
  out <- patch("--install-from", "_local_data/cube_patch_all", "--dry-run")
  expect_null(attr(out, "status"))
  expect_equal(tools::md5sum(file.path(w$cdir, c("site_1_One_data.csv", "site_2_Two_data.csv"))), orig)   # still untouched
  out <- patch("--install-from", "_local_data/cube_patch_all")
  expect_null(attr(out, "status"))
  bk <- file.path(w$dir, "_local_data", paste0("cube_backup_", format(Sys.Date(), "%Y%m%d")))
  expect_equal(unname(tools::md5sum(file.path(bk, "ramsar_site_data_100m_europe/Belgium", c("site_1_One_data.csv", "site_2_Two_data.csv")))), unname(orig))
  expect_true(file.exists(file.path(bk, "ramsar_site_data_100m_unknown/x/site_1_One_data.csv")))
  expect_false(file.exists(file.path(unk, "site_1_One_data.csv")))
  il <- list.files(file.path(w$dir, "_local_data/cube_patch_work"), pattern = "^install_log_[0-9_]+\\.csv$", full.names = TRUE)
  expect_length(il, 1)
  il <- fread(il)
  expect_true(all(il$verified)); expect_setequal(il$action, "replaced")
  rd <- fread(file.path(w$dir, "_local_data/cube_patch_work/update_files_read.csv"))
  expect_equal(rd$rows_other_years, c(2, 0)); expect_equal(rd$rows_other_taxa, c(0, 1))
  run_check(w, "--step", "all", "--tag", "patched")
  st <- fread(file.path(w$dir, "check/compare_patched/sites.csv"))[site_id %in% c("site_1", "site_2")]
  expect_equal(st$rows_cube, c(nrow(w$site1), nrow(w$site2)))
  expect_equal(st$rows_same, st$rows_ref); expect_equal(st$rows_same, st$rows_cube)
})

test_that("rows that differ only in the higher-taxon counts are counted apart", {
  row <- function(sp, occ, kc) data.table(kingdom = "Animalia", kingdomkey = "1", phylum = "Chordata", phylumkey = "44", class = "Aves",
                                          classkey = "212", order = "O", orderkey = "1", family = "F", familykey = "1", genus = "G",
                                          genuskey = "1", species = sp, specieskey = sp, year = "2010", mgrscellcode = "31UES000000",
                                          kingdomcount = kc, phylumcount = kc, classcount = "1", ordercount = "1", familycount = "1",
                                          genuscount = "1", occurrences = occ, mintemporaluncertainty = "60",
                                          mincoordinateuncertaintyinmeters = "100")
  ref <- chk_normalise(rbind(row("A", "1", "5"), row("B", "2", "5"), row("C", "1", "5")))
  cube <- chk_normalise(rbind(row("A", "1", "1"), row("B", "2", "5"), row("C", "3", "5")))
  ll <- mgrs::mgrs_to_latlng("31UES000000")
  box <- sf::st_as_sfc(sf::st_bbox(c(xmin = ll$lng - 0.01, ymin = ll$lat - 0.01, xmax = ll$lng + 0.01, ymax = ll$lat + 0.01), crs = 4326))
  r <- chk_compare_site(cube, ref, CHK_COLS, sf::st_sf(site_id = "site_1", geometry = box), 2024)$counts
  expect_equal(r$rows_same, 2)            # A (counts only) and B
  expect_equal(r$count_cols_differ, 1)    # A
  expect_equal(r$missing_before, 1)       # C with 1 occurrence
  expect_equal(r$extra_inside_site, 1)    # C with 3 occurrences
})

test_that("older rows rebuilt from the reference pass the check", {
  skip_if_not_installed("withr")
  w <- world()
  repo <- normalizePath(file.path(here, ".."))
  run_check(w, "--step", "inventory"); run_check(w, "--step", "reference")
  file.symlink(file.path(repo, "scripts"), file.path(w$dir, "scripts"))
  rs <- function(script, ...) withr::with_dir(w$dir, system2("Rscript", c(script, ...), stdout = TRUE, stderr = TRUE))
  # site 2 filed under the wrong continent (its rows come from the Europe cube): moved back by the rebuild
  adir <- file.path(w$dir, "inst/extdata/ramsar_site_data_100m_asia/Belgium"); dir.create(adir, recursive = TRUE)
  file.rename(file.path(w$cdir, "site_2_Two_data.csv"), file.path(adir, "site_2_Two_data.csv"))
  orig <- tools::md5sum(c(file.path(w$cdir, "site_1_One_data.csv"), file.path(adir, "site_2_Two_data.csv")))
  out <- rs("scripts/rebuild_site_cubes_from_reference.R", "--check-dir", "check", "--min-rows", "1", "--cores", "1", "--dry-run")
  expect_null(attr(out, "status"))
  expect_false(dir.exists(file.path(w$dir, "_local_data/cube_rebuild")))
  out <- rs("scripts/rebuild_site_cubes_from_reference.R", "--check-dir", "check", "--min-rows", "1", "--cores", "1")
  expect_null(attr(out, "status"))
  expect_equal(tools::md5sum(c(file.path(w$cdir, "site_1_One_data.csv"), file.path(adir, "site_2_Two_data.csv"))), orig)   # cubes untouched
  sm <- fread(file.path(w$dir, "_local_data/cube_rebuild_work/rebuild_summary.csv"))
  expect_equal(sm[site_id == "site_1", .(action, missing_before, extra_before)], data.table(action = "rebuilt", missing_before = 3L, extra_before = 2L))
  expect_equal(sm[site_id == "site_2", .(action, cube, moved_from)],
               data.table(action = "moved", cube = "ramsar_site_data_100m_europe/Belgium/site_2_Two_data.csv",
                          moved_from = "ramsar_site_data_100m_asia/Belgium/site_2_Two_data.csv"))
  expect_equal(sm[site_id == "site_3", .(action, cube)],
               data.table(action = "new cube", cube = "ramsar_site_data_100m_europe/Belgium/site_3_Three_data.csv"))
  # rerunning with other settings in the same work folder is refused
  out <- suppressWarnings(rs("scripts/rebuild_site_cubes_from_reference.R", "--check-dir", "check", "--min-rows", "5", "--cores", "1"))
  expect_false(is.null(attr(out, "status")))
  out <- rs("scripts/patch_site_cubes_update.R", "--install-from", "_local_data/cube_rebuild", "--work-dir", "_local_data/cube_rebuild_work")
  expect_null(attr(out, "status"))
  run_check(w, "--step", "compare", "--tag", "rebuilt")
  st <- fread(file.path(w$dir, "check/compare_rebuilt/sites.csv"))
  s1 <- st[site_id == "site_1"]
  expect_equal(c(s1$missing_before, s1$extra_inside_site, s1$extra_outside_site, s1$extra_before), c(0, 0, 0, 0))
  expect_equal(s1$missing_from, 1)                     # the 2025 snake: rows from 2024 on are kept as they were
  s3 <- st[site_id == "site_3"]
  expect_equal(c(s3$n_cube_files, s3$rows_cube, s3$rows_same), c(1, 1, 1))
  expect_equal(st[site_id == "site_2", .(n_cube_files, rows_same)], data.table(n_cube_files = 1L, rows_same = nrow(w$site2)))
  expect_false(file.exists(file.path(adir, "site_2_Two_data.csv")))      # the misfiled copy is in the backup folder
  bk <- file.path(w$dir, "_local_data", paste0("cube_backup_", format(Sys.Date(), "%Y%m%d")), "ramsar_site_data_100m_asia/Belgium/site_2_Two_data.csv")
  expect_equal(unname(tools::md5sum(bk)), unname(orig[2]))
  # --no-refile leaves misfiled cubes where they are
  file.copy(bk, file.path(adir, "site_2_Two_data.csv")); file.remove(file.path(w$cdir, "site_2_Two_data.csv"))
  out <- rs("scripts/rebuild_site_cubes_from_reference.R", "--check-dir", "check", "--min-rows", "1", "--cores", "1", "--no-refile",
            "--out-dir", "_local_data/cube_rebuild2", "--work-dir", "_local_data/cube_rebuild2_work")
  expect_null(attr(out, "status"))
  sm2 <- fread(file.path(w$dir, "_local_data/cube_rebuild2_work/rebuild_summary.csv"))
  expect_equal(sm2[site_id == "site_2", action], "kept")
})
