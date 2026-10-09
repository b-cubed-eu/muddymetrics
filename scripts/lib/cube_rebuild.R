# Rebuilding the older rows of site cubes from the rows the cube check has
# collected (scripts/rebuild_site_cubes_from_reference.R).
#
# The reference step of scripts/check_site_cubes.R reads every line of the
# continental cubes and the update files and keeps, per site, every row whose
# cell lies in the site's outline in the Ramsar shapefile. Where a site cube
# lacks many of these rows (cubes split incompletely in June 2026, outlines that
# are invalid in planar geometry, rows filed under a neighbouring continent),
# its rows before `from_year` are replaced by the reference rows. Rows from
# `from_year` on are kept as they are in the cube (they were rebuilt by
# scripts/patch_site_cubes_update.R and agree with the check).
#
# Requires data.table and the functions of scripts/lib/cube_check.R.

#' Reference rows of one bucket, with the source each row came from (the name
#' of the continental cube or update file, without the part number).
rb_read_reference_bucket <- function(ref_dir, bucket) {
  dirs <- list.dirs(ref_dir, recursive = FALSE)
  dirs <- dirs[file.exists(file.path(dirs, "done.csv"))]
  fs <- file.path(dirs, sprintf("bucket_%03d.tsv", bucket))
  keep <- file.exists(fs)
  if (!any(keep)) return(data.table::data.table(matrix(character(), 0, length(CHK_COLS) + 2,
                                                       dimnames = list(NULL, c("site_id", CHK_COLS, "source")))))
  data.table::rbindlist(Map(function(f, d) {
    x <- data.table::fread(f, sep = "\t", quote = "", colClasses = "character", encoding = "UTF-8", showProgress = FALSE)
    x[, source := sub("_part[0-9]+$", "", basename(d))]
    x
  }, fs[keep], dirs[keep]))
}

#' Keys for matching rows: every cube column except the higher-taxon counts,
#' normalised as in the check.
rb_keys <- function(dt) {
  n <- chk_normalise(dt, CHK_COLS)
  do.call(paste, c(as.list(n[, setdiff(CHK_COLS, CHK_COUNT_COLS), with = FALSE]), sep = CHK_SEP))
}

#' How many of the reference rows before `from_year` are missing from the cube
#' (matched on rb_keys(), as a multiset), and how many cube rows before
#' `from_year` are not in the reference.
rb_missing_before <- function(cube, ref, from_year) {
  yc <- suppressWarnings(as.numeric(cube$year)); yr <- suppressWarnings(as.numeric(ref$year))
  cb <- cube[!is.na(yc) & yc < from_year]; rf <- ref[!is.na(yr) & yr < from_year]
  tc <- data.table::data.table(k = rb_keys(cb))[, .(n_cube = .N), by = k]
  tr <- data.table::data.table(k = rb_keys(rf))[, .(n_ref = .N), by = k]
  m <- merge(tc, tr, by = "k", all = TRUE)
  m[is.na(n_cube), n_cube := 0L]; m[is.na(n_ref), n_ref := 0L]
  list(rows_cube = nrow(cb), rows_ref = nrow(rf), missing = sum(pmax(m$n_ref - m$n_cube, 0L)),
       extra = sum(pmax(m$n_cube - m$n_ref, 0L)))
}

#' New contents of a site cube: the reference rows before `from_year`
#' followed by the cube's own rows from `from_year` on (for a site without a
#' cube, the reference rows from `from_year` on). Cube columns only, as text.
rb_rebuilt_rows <- function(cube, ref, from_year) {
  yr <- suppressWarnings(as.numeric(ref$year))
  before <- ref[!is.na(yr) & yr < from_year, CHK_COLS, with = FALSE]
  from <- if (is.null(cube)) ref[!is.na(yr) & yr >= from_year, CHK_COLS, with = FALSE] else {
    yc <- suppressWarnings(as.numeric(cube$year))
    cube[!is.na(yc) & yc >= from_year, CHK_COLS, with = FALSE]
  }
  rbind(before, from)
}

#' Continent folder a site's cube belongs in: the continent whose continental
#' cube holds most of the site's reference rows (years before 2024). NA when
#' there are no such rows or two continents hold the same number.
rb_continent_of <- function(sources) {
  s <- sources[grepl("_100m_mgrs_all$", sources)]
  if (!length(s)) return(NA_character_)
  tb <- sort(table(sub("_100m_mgrs_all$", "", s)), decreasing = TRUE)
  if (length(tb) > 1 && tb[[1]] == tb[[2]]) return(NA_character_)
  names(tb)[1]
}
