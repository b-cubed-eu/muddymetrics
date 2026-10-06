# Functions for patching site cubes with a newer GBIF occurrence cube.
#
# Used by scripts/patch_site_cubes_update.R and its tests. A site cube is the
# per-site CSV under inst/extdata/ramsar_site_data_100m_<continent>/<country>/;
# an update cube is a GBIF SQL download with the same 25 columns, one row per
# species x 100 m MGRS cell x year.
#
# Cells are assigned to sites the same way as in the original cube split
# (scripts/split_data_cubes_*.R): the MGRS cell code is converted to a point
# with mgrs::mgrs_to_latlng() and the point is intersected with the site
# outline (planar geometry, sf_use_s2(FALSE)).

CUBE_COLS <- c("kingdom", "kingdomkey", "phylum", "phylumkey", "class", "classkey",
               "order", "orderkey", "family", "familykey", "genus", "genuskey",
               "species", "specieskey", "year", "mgrscellcode", "kingdomcount",
               "phylumcount", "classcount", "ordercount", "familycount", "genuscount",
               "occurrences", "mintemporaluncertainty", "mincoordinateuncertaintyinmeters")

CONTINENTS <- c("africa", "antarctica", "asia", "europe", "northamerica", "oceania", "southamerica")

MGRS_100M <- "^[0-9]{1,2}[C-X][A-Z]{2}[0-9]{6}$"

#' Site id ("site_123") from a cube or outline file name. Unlike splitting the
#' name on ".", this is safe for site names that contain dots ("St.").
site_id_of <- function(path) {
  b <- basename(path)
  ifelse(grepl("^site_[0-9]+", b), sub("^(site_[0-9]+).*$", "\\1", b), NA_character_)
}

#' Every site cube file under the continent folders (all files, not one per
#' site), with continent, country and site id taken from the path.
list_site_cubes <- function(root = "inst/extdata", continents = CONTINENTS) {
  out <- list()
  for (cont in continents) {
    d <- file.path(root, paste0("ramsar_site_data_100m_", cont))
    if (!dir.exists(d)) next
    fs <- list.files(d, pattern = "_data\\.csv$", recursive = TRUE, full.names = TRUE)
    if (!length(fs)) next
    out[[cont]] <- data.frame(continent = cont, country = basename(dirname(fs)),
                              site_id = vapply(fs, site_id_of, ""), file_path = fs,
                              stringsAsFactors = FALSE, row.names = NULL)
  }
  do.call(rbind, out)
}

#' Outline file for a site id: the .wkt whose name starts with "<site_id>_",
#' preferring the one in `country`'s folder when there is more than one.
find_site_wkt <- function(site_id, wkt_files, country = NULL) {
  hit <- wkt_files[startsWith(basename(wkt_files), paste0(site_id, "_"))]
  if (length(hit) > 1 && !is.null(country)) {
    same <- hit[basename(dirname(hit)) == country]
    if (length(same)) hit <- same
  }
  if (length(hit)) hit[1] else NA_character_
}

#' Read site outlines into an sf object with one row per site.
read_site_outlines <- function(site_ids, wkt_paths) {
  geoms <- lapply(wkt_paths, function(p) sf::st_as_sfc(paste(readLines(p, warn = FALSE), collapse = ""), crs = 4326)[[1]])
  s <- sf::st_sf(site_id = site_ids, geometry = sf::st_sfc(geoms, crs = 4326))
  suppressWarnings(sf::st_make_valid(s))
}

#' 100 km MGRS squares ("30UVB") touched by a site outline: a 21 x 21 grid of
#' points over the bounding box plus all vertices, as in the original split.
site_squares <- function(geom) {
  bb <- sf::st_bbox(geom)
  xs <- seq(bb[["xmin"]], bb[["xmax"]], by = max(0.01, (bb[["xmax"]] - bb[["xmin"]]) / 20))
  ys <- seq(bb[["ymin"]], bb[["ymax"]], by = max(0.01, (bb[["ymax"]] - bb[["ymin"]]) / 20))
  pts <- rbind(expand.grid(lng = xs, lat = ys),
               data.frame(lng = sf::st_coordinates(geom)[, "X"], lat = sf::st_coordinates(geom)[, "Y"]))
  pts <- unique(round(pts, 4))
  sq <- tryCatch(mapply(mgrs::latlng_to_mgrs, pts$lat, pts$lng, MoreArgs = list(precision = 0)),
                 error = function(e) character())
  unique(stats::na.omit(unlist(sq)))
}

#' Sites (site_id) whose outline contains each MGRS 100 m cell. Returns a
#' data.table(mgrscellcode, site_id); a cell can belong to several sites.
assign_cells <- function(cells, outlines) {
  cells <- unique(cells[grepl(MGRS_100M, cells)])
  if (!length(cells)) return(data.table::data.table(mgrscellcode = character(), site_id = character()))
  ll <- mgrs::mgrs_to_latlng(cells)
  pts <- sf::st_as_sf(data.frame(mgrscellcode = ll$mgrs, lng = ll$lng, lat = ll$lat), coords = c("lng", "lat"), crs = 4326)
  s2 <- sf::sf_use_s2(); on.exit(suppressMessages(sf::sf_use_s2(s2)))
  suppressMessages(sf::sf_use_s2(FALSE))
  j <- suppressMessages(sf::st_join(pts, outlines, join = sf::st_intersects, left = FALSE))
  data.table::data.table(mgrscellcode = j$mgrscellcode, site_id = j$site_id)
}

#' Does a site cube still need the update? TRUE unless it has rows from
#' `from_year` on and most of them carry new-backbone species keys. Old-backbone
#' keys are numeric; new-backbone keys are short codes such as "7CYVX", a few of
#' which happen to be all digits, hence "most" rather than "all".
#' Files with no rows from `from_year` on also return TRUE (nothing to lose).
needs_patch <- function(path, from_year = 2024) {
  needs_patch_rows(data.table::fread(path, select = c("year", "specieskey"), colClasses = "character", showProgress = FALSE), from_year)
}
needs_patch_rows <- function(x, from_year = 2024) {
  y <- suppressWarnings(as.integer(x$year))
  k <- x$specieskey[!is.na(y) & y >= from_year]
  !length(k) || mean(grepl("^[0-9]+$", k)) > 0.5
}

#' Stream an update cube in chunks of lines, keeping rows whose cell lies in one
#' of `squares` (5-character 100 km squares). Calls on_chunk(dt) for each chunk
#' with matching rows (all columns character). Returns the number of data lines read.
stream_update <- function(path, squares, on_chunk, chunk_lines = 2e6, progress = TRUE) {
  # Read raw bytes: a connection opened with encoding = "UTF-8" stops silently at
  # the first invalid byte sequence. Invalid sequences are dropped line by line.
  con <- file(path, open = "rb")
  on.exit(close(con))
  first <- readLines(con, n = 1, warn = FALSE)
  header <- strsplit(sub("\r$", "", first), "\t", fixed = FALSE)[[1]]
  if (!identical(header, CUBE_COLS)) stop("Unexpected columns in ", path, ": ", paste(header, collapse = ","))
  n <- 0; bytes <- nchar(first, type = "bytes") + 1; t0 <- Sys.time()
  repeat {
    lines <- readLines(con, n = chunk_lines, warn = FALSE, skipNul = TRUE)
    if (!length(lines)) break
    n <- n + length(lines)
    bytes <- bytes + sum(nchar(lines, type = "bytes", allowNA = TRUE), na.rm = TRUE) + length(lines)
    lines <- iconv(lines, "UTF-8", "UTF-8", sub = "")
    lines[is.na(lines)] <- ""
    # cheap pre-filter on the cell column (16th field) before parsing
    cell <- stringi::stri_match_first_regex(lines, "^(?:[^\t]*\t){15}([^\t]*)\t")[, 2]
    keep <- !is.na(cell) & substr(cell, 1, 5) %in% squares &
      stringi::stri_count_fixed(lines, "\t") == length(CUBE_COLS) - 1   # whole lines only
    if (any(keep)) {
      dt <- data.table::fread(text = paste0(paste(lines[keep], collapse = "\n"), "\n"), sep = "\t", header = FALSE, quote = "",
                              colClasses = "character", col.names = CUBE_COLS, fill = TRUE, showProgress = FALSE)
      dt <- dt[!is.na(dt$mincoordinateuncertaintyinmeters) & nzchar(dt$mincoordinateuncertaintyinmeters)]
      on_chunk(dt)
    }
    if (progress) cat(sprintf("  %s lines read (%.1f min)\n", format(n, big.mark = ",", scientific = FALSE), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  size <- file.size(path)
  if (bytes < 0.999 * size) stop(sprintf("Stopped after %s of %s bytes of %s: the update was not read to the end.",
                                         format(bytes, big.mark = ","), format(size, big.mark = ","), path))
  n
}

#' New contents of a site cube: its rows before `from_year` (cube columns only,
#' extra site columns dropped) followed by the update rows for the site.
patched_rows <- function(old, new, from_year = 2024) {
  old <- data.table::as.data.table(old)
  miss <- setdiff(CUBE_COLS, names(old))
  if (length(miss)) stop("Site cube lacks columns: ", paste(miss, collapse = ", "))
  keep <- old[suppressWarnings(as.integer(old$year)) < from_year, CUBE_COLS, with = FALSE]
  new <- data.table::as.data.table(new)[, CUBE_COLS, with = FALSE]
  for (cc in CUBE_COLS) { data.table::set(keep, j = cc, value = as.character(keep[[cc]])); data.table::set(new, j = cc, value = as.character(new[[cc]])) }
  rbind(keep, new)
}

#' Rows from `from_year` on as comparable keys (cell, year, species key, occurrences).
row_keys <- function(dt, from_year = 2024) {
  dt <- data.table::as.data.table(dt)
  y <- suppressWarnings(as.integer(dt$year))
  d <- dt[!is.na(y) & y >= from_year]
  sort(paste(d$mgrscellcode, d$year, d$specieskey, as.numeric(d$occurrences), sep = "|"))
}
