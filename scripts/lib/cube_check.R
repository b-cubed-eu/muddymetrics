# Independent checks of the site cubes against their sources.
#
# Used by scripts/check_site_cubes.R and its tests. The checks rebuild, from the
# source files alone, what every site cube should contain, and compare that
# row by row with the site cubes on disk. They are written separately from the
# code that built the cubes (scripts/split_data_cubes_*.R,
# scripts/patch_site_cubes_update.R, scripts/lib/cube_update.R), so that an
# error in that code is not repeated here:
#
#   - site outlines come from the Ramsar shapefile, not from the .wkt files;
#   - source files are read in byte ranges with their own reader, and every
#     line is accounted for (lines that are not 25 tab-separated fields are
#     counted, not dropped silently);
#   - cells are not pre-selected by 100 km square: every valid cell is placed
#     on the map and tested against every outline.
#
# The rule for which cells belong to a site is the one the cubes were built
# with: a 100 m MGRS cell belongs to a site when the point returned by
# mgrs::mgrs_to_latlng() for it (the south-west corner of the cell) lies in
# or on the site outline, with planar geometry (sf_use_s2(FALSE)).
#
# Requires data.table, sf, mgrs and stringi.

CHK_COLS <- c("kingdom", "kingdomkey", "phylum", "phylumkey", "class", "classkey",
              "order", "orderkey", "family", "familykey", "genus", "genuskey",
              "species", "specieskey", "year", "mgrscellcode", "kingdomcount",
              "phylumcount", "classcount", "ordercount", "familycount", "genuscount",
              "occurrences", "mintemporaluncertainty", "mincoordinateuncertaintyinmeters")
CHK_NUM_COLS <- c("year", "kingdomcount", "phylumcount", "classcount", "ordercount", "familycount",
                  "genuscount", "occurrences", "mintemporaluncertainty", "mincoordinateuncertaintyinmeters")
CHK_TAXON_COLS <- c("kingdom", "kingdomkey", "phylum", "phylumkey", "class", "classkey",
                    "order", "orderkey", "family", "familykey", "genus", "genuskey",
                    "species", "specieskey")
# Columns compared once the taxonomy has been harmonised: the per-rank count
# columns are not kept meaningfully when rows are merged, so they are left out.
CHK_HARMONISED_COLS <- c(CHK_TAXON_COLS, "year", "mgrscellcode", "occurrences",
                         "mintemporaluncertainty", "mincoordinateuncertaintyinmeters")
# Occurrences of the row's kingdom, phylum, ... genus in the cell and year. They
# are not used by the screen or the indicators, and in the older site cubes
# they are lower than in the continental cubes; rows that differ only in these
# columns are counted separately (count_cols_differ), not as differing rows.
CHK_COUNT_COLS <- c("kingdomcount", "phylumcount", "classcount", "ordercount", "familycount", "genuscount")
# 100 m MGRS cell: zone, latitude band (C-X without I and O), 100 km square
# letters (without I and O), 3 + 3 digits.
CHK_MGRS_100M <- "^[0-9]{1,2}[C-HJ-NP-X][A-HJ-NP-Z]{2}[0-9]{6}$"
CHK_SEP <- "\x1f"

chk_site_id_of <- function(path) {
  b <- basename(path)
  ifelse(grepl("^site_[0-9]+_", b), sub("^(site_[0-9]+)_.*$", "\\1", b), NA_character_)
}

# ---------------------------------------------------------------------------
# Outlines and inventory
# ---------------------------------------------------------------------------

#' Site outlines from the Ramsar shapefile: one row per ramsarid (parts
#' unioned after repair with st_make_valid), as sf with site_id, name,
#' country, iso3, official area (ha).
chk_read_outlines <- function(shp) {
  s <- suppressWarnings(sf::st_read(shp, quiet = TRUE))
  old <- sf::sf_use_s2(); on.exit(suppressMessages(sf::sf_use_s2(old)))
  suppressMessages(sf::sf_use_s2(FALSE))
  s <- suppressWarnings(sf::st_make_valid(s))
  ids <- as.character(s$ramsarid)
  idx <- split(seq_len(nrow(s)), ids)
  to_utf8 <- function(x) { x <- as.character(x); bad <- !stringi::stri_enc_isutf8(x); x[bad] <- iconv(x[bad], "latin1", "UTF-8"); x }
  geoms <- lapply(idx, function(i) {
    g <- suppressMessages(sf::st_union(sf::st_geometry(s)[i]))
    if (!all(sf::st_geometry_type(g) %in% c("POLYGON", "MULTIPOLYGON")))
      g <- sf::st_collection_extract(g, "POLYGON")
    sf::st_cast(g, "MULTIPOLYGON")[[1]]
  })
  first <- vapply(idx, `[`, 1L, 1)
  sf::st_sf(site_id = paste0("site_", names(idx)), name = to_utf8(s$officialna[first]),
            country = to_utf8(s$country_en[first]), iso3 = as.character(s$iso3[first]),
            area_off_ha = as.numeric(s$area_off[first]), n_parts = unname(lengths(idx)),
            geometry = sf::st_sfc(unname(geoms), crs = sf::st_crs(s)))
}

#' Site cube files: every *_data.csv under the ramsar_site_data_100m_* folders
#' of `root`, except the folders holding the 2024-2026 update.
chk_list_cubes <- function(root) {
  tops <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  tops <- tops[grepl("^ramsar_site_data_100m_", basename(tops)) & !grepl("world", basename(tops))]
  out <- data.table::rbindlist(lapply(tops, function(top) {
    fs <- list.files(top, pattern = "_data\\.csv$", recursive = TRUE, full.names = TRUE)
    data.table::data.table(file = fs, folder = rep(sub("^ramsar_site_data_100m_", "", basename(top)), length(fs)),
                           country_dir = basename(dirname(fs)), site_id = chk_site_id_of(fs))
  }))
  if (!nrow(out)) out <- data.table::data.table(file = character(), folder = character(), country_dir = character(), site_id = character())
  out
}

#' One row per site in the shapefile, and one per cube or outline file whose
#' site is not in the shapefile, with what was found for it.
chk_inventory <- function(outlines, cubes, wkt_files) {
  o <- data.table::data.table(site_id = outlines$site_id, name = outlines$name, shp_country = outlines$country)
  w <- data.table::data.table(wkt = wkt_files, site_id = chk_site_id_of(wkt_files), wkt_country = basename(dirname(wkt_files)))
  cn <- cubes[, .(n_cubes = .N, cube_files = paste(file, collapse = " | "), cube_folders = paste(unique(folder), collapse = ","),
                  cube_country = paste(unique(country_dir), collapse = ",")), by = site_id]
  wn <- w[, .(n_wkt = .N, wkt_files = paste(wkt, collapse = " | "), wkt_country = paste(unique(wkt_country), collapse = ",")), by = site_id]
  inv <- merge(merge(o, cn, by = "site_id", all = TRUE), wn, by = "site_id", all = TRUE)
  inv[is.na(n_cubes), n_cubes := 0L]; inv[is.na(n_wkt), n_wkt := 0L]
  inv[, in_shapefile := site_id %in% outlines$site_id]
  inv[, `:=`(issues = "", notes = "")]
  add <- function(cond, txt) inv[cond, issues := ifelse(issues == "", txt, paste(issues, txt, sep = "; "))]
  note <- function(cond, txt) inv[cond, notes := ifelse(notes == "", txt, paste(notes, txt, sep = "; "))]
  add(!inv$in_shapefile, "not in the shapefile")
  add(inv$in_shapefile & inv$n_cubes == 0, "no cube file")
  add(inv$n_cubes > 1, "more than one cube file")
  add(inv$in_shapefile & inv$n_wkt == 0, "no outline file")
  add(inv$n_wkt > 1, "more than one outline file")
  add(grepl("unknown", inv$cube_folders), "cube in the unknown folder")
  # Folder names are the country names of the outline export (e.g. "Turkey" where the shapefile now has
  # "Türkiye"): a different name is noted, not counted as an issue.
  note(inv$n_cubes > 0 & !is.na(inv$shp_country) & inv$cube_country != inv$shp_country, "cube folder name differs from the shapefile's country name")
  note(inv$n_wkt > 0 & !is.na(inv$shp_country) & inv$wkt_country != inv$shp_country, "outline folder name differs from the shapefile's country name")
  note(inv$n_cubes > 0 & inv$n_wkt > 0 & inv$cube_country != inv$wkt_country, "cube and outline in folders of different countries")
  inv[order(as.integer(sub("site_", "", site_id)))]
}

#' Compare each .wkt outline (as the cube split read it: st_make_valid of the
#' text) with the shapefile outline. rel_symdiff is the area that lies in one
#' outline but not the other, as a share of the shapefile outline (planar
#' areas in degrees squared). mean_shift_m is that area divided by the length
#' of the shapefile boundary, converted to metres at the site's latitude: the
#' average distance between the two boundaries (a few metres for rounded
#' coordinates, much more for a different outline).
chk_compare_wkt <- function(outlines, wkt_files) {
  old <- sf::sf_use_s2(); on.exit(suppressMessages(sf::sf_use_s2(old)))
  suppressMessages(sf::sf_use_s2(FALSE))
  ids <- chk_site_id_of(wkt_files)
  rows <- lapply(seq_along(wkt_files), function(i) {
    j <- match(ids[i], outlines$site_id)
    base <- data.table::data.table(site_id = ids[i], wkt = wkt_files[i])
    if (is.na(j)) return(cbind(base, status = "not in the shapefile"))
    res <- tryCatch({
      g <- sf::st_as_sfc(paste(readLines(wkt_files[i], warn = FALSE), collapse = ""), crs = sf::st_crs(outlines))
      g <- suppressWarnings(sf::st_make_valid(g))
      if (!all(sf::st_geometry_type(g) %in% c("POLYGON", "MULTIPOLYGON"))) g <- sf::st_collection_extract(g, "POLYGON")
      s <- sf::st_geometry(outlines)[j]
      planar <- function(x) sf::st_set_crs(x, NA)                                  # no extra packages needed
      parea <- function(x) as.numeric(sum(sf::st_area(planar(x))))
      a_s <- parea(s); a_w <- parea(g)
      sd <- parea(suppressMessages(sf::st_sym_difference(sf::st_union(g), s)))
      perim <- as.numeric(sum(sf::st_length(sf::st_cast(planar(s), "MULTILINESTRING"))))
      bw <- sf::st_bbox(g); bs <- sf::st_bbox(s)
      lat <- (bs[["ymin"]] + bs[["ymax"]]) / 2
      km_per_deg <- 111.32 * sqrt(cos(lat * pi / 180))                            # geometric mean of the two axes
      data.table::data.table(status = "compared", area_ratio = a_w / a_s, rel_symdiff = sd / a_s,
                             mean_shift_m = 1000 * km_per_deg * sd / perim,
                             bbox_shift_deg = max(abs(as.numeric(bw) - as.numeric(bs))),
                             n_vertices_wkt = nrow(sf::st_coordinates(g)), n_vertices_shp = nrow(sf::st_coordinates(s)))
    }, error = function(e) data.table::data.table(status = paste("could not compare:", conditionMessage(e))))
    cbind(base, res)
  })
  data.table::rbindlist(rows, fill = TRUE)
}

# ---------------------------------------------------------------------------
# Reading source files in byte ranges
# ---------------------------------------------------------------------------

#' Bytes per line ending (1 for LF, 2 for CRLF), from the first line.
chk_eol_bytes <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  b <- readBin(con, "raw", n = 1e6)
  nl <- which(b == as.raw(10))[1]
  if (is.na(nl)) stop("No line end in the first 1 MB of ", path)
  if (nl > 1 && b[nl - 1] == as.raw(13)) 2L else 1L
}

#' Does the file end with a line end?
chk_final_eol <- function(path) {
  size <- file.size(path)
  con <- file(path, "rb"); on.exit(close(con))
  seek(con, size - 1)
  identical(readBin(con, "raw", n = 1), as.raw(10))
}

#' Split a file into byte ranges of about `part_bytes`.
chk_plan_parts <- function(path, part_bytes) {
  size <- file.size(path)
  part_bytes <- max(1, floor(part_bytes))       # whole bytes
  starts <- seq(0, size - 1, by = part_bytes)
  data.table::data.table(path = path, part = seq_along(starts), start = starts, end = pmin(starts + part_bytes, size), size = size)
}

#' Read the lines that start in [start, end) of a file and pass them to
#' on_chunk() in chunks. A line belongs to the range that contains its first
#' byte, so consecutive ranges together read every line exactly once. The
#' header (first line) is checked against CHK_COLS and not passed on.
#' Returns the byte position of the first line read, the position after the
#' last one, and the number of lines.
chk_read_range <- function(path, start, end, on_chunk, chunk_lines = 1e6) {
  eol <- chk_eol_bytes(path)
  con <- file(path, open = "rb")
  on.exit(close(con))
  if (start == 0) {
    h <- readLines(con, n = 1, warn = FALSE)
    hdr <- strsplit(sub("\r$", "", h), "\t", fixed = TRUE)[[1]]
    if (!identical(hdr, CHK_COLS)) stop("Unexpected header in ", path, ": ", paste(hdr, collapse = ","))
    pos <- nchar(h, type = "bytes") + eol
  } else {
    seek(con, where = start - 1)
    rest <- readLines(con, n = 1, warn = FALSE)        # end of the line that holds byte start - 1
    pos <- start - 1 + nchar(rest, type = "bytes") + eol
  }
  first <- pos; n <- 0
  while (pos < end) {
    lines <- readLines(con, n = chunk_lines, warn = FALSE)
    if (!length(lines)) break
    len <- nchar(lines, type = "bytes", allowNA = TRUE) + eol
    starts <- pos + c(0, cumsum(len)[-length(len)])
    own <- starts < end
    lines <- lines[own]
    pos <- pos + sum(len[own])
    n <- n + length(lines)
    on_chunk(lines)
    if (!all(own)) break
  }
  list(first = first, after = pos, lines = n)
}

#' Check that the ranges of a file were read without gaps or overlaps and up to
#' the end of the file. `parts`: data.table with part, first, after (from
#' chk_read_range), size, eol_bytes and final_eol.
chk_check_tiling <- function(parts) {
  p <- parts[order(part)]
  expected_end <- p$size[1] + if (isTRUE(p$final_eol[1])) 0 else p$eol_bytes[1]
  gaps <- if (nrow(p) > 1) which(p$first[-1] != p$after[-nrow(p)]) else integer()
  list(ok = !length(gaps) && p$after[nrow(p)] == expected_end,
       gaps_after_part = p$part[gaps], end_read = p$after[nrow(p)], end_expected = expected_end)
}

#' Split raw lines into a character data.table with the 25 cube columns.
#' Returns the table plus counts of lines that were empty, did not have 25
#' fields, or had invalid UTF-8 (those bytes are removed, as in the cubes).
chk_parse_lines <- function(lines) {
  empty <- lines == ""
  nf <- stringi::stri_count_fixed(lines, "\t") + 1L
  ok <- !empty & nf == length(CHK_COLS)
  bad_examples <- utils::head(lines[!ok & !empty], 3)
  lines <- lines[ok]
  utf <- stringi::stri_enc_isutf8(lines)
  if (any(!utf)) lines[!utf] <- iconv(lines[!utf], "UTF-8", "UTF-8", sub = "")
  Encoding(lines) <- "UTF-8"          # the bytes are UTF-8 whatever the session's locale
  m <- stringi::stri_split_fixed(lines, "\t", n = length(CHK_COLS), simplify = TRUE)
  if (!length(lines)) m <- matrix(character(), 0, length(CHK_COLS))
  dt <- data.table::as.data.table(m)
  data.table::setnames(dt, CHK_COLS)
  list(rows = dt, n_empty = sum(empty), n_wrong_fields = sum(!ok & !empty), n_bad_utf8 = sum(!utf),
       bad_examples = bad_examples)
}

# ---------------------------------------------------------------------------
# Cells to sites
# ---------------------------------------------------------------------------

#' Sites whose outline contains the point of each cell. Returns
#' data.table(mgrscellcode, site_id) (a cell can belong to several sites) and
#' the counts of cells with an invalid code or no coordinates.
chk_cells_to_sites <- function(cells, outlines) {
  cells <- unique(cells)
  valid <- grepl(CHK_MGRS_100M, cells)
  ll <- if (any(valid)) suppressWarnings(mgrs::mgrs_to_latlng(cells[valid])) else data.frame(mgrs = character(), lat = numeric(), lng = numeric())
  if (!identical(as.character(ll$mgrs), cells[valid])) stop("mgrs_to_latlng() did not return the cells in input order")
  conv <- !is.na(ll$lat) & !is.na(ll$lng)
  out <- data.table::data.table(mgrscellcode = character(), site_id = character())
  if (any(conv)) {
    pts <- sf::st_as_sf(data.frame(mgrscellcode = cells[valid][conv], lng = ll$lng[conv], lat = ll$lat[conv]),
                        coords = c("lng", "lat"), crs = sf::st_crs(outlines))
    old <- sf::sf_use_s2(); on.exit(suppressMessages(sf::sf_use_s2(old)))
    suppressMessages(sf::sf_use_s2(FALSE))
    hit <- suppressMessages(sf::st_intersects(pts, outlines))
    n <- lengths(hit)
    if (sum(n)) out <- data.table::data.table(mgrscellcode = rep(pts$mgrscellcode, n), site_id = outlines$site_id[unlist(hit)])
  }
  list(pairs = out, n_invalid_code = sum(!valid), n_no_coords = sum(!conv))
}

# ---------------------------------------------------------------------------
# Comparing a cube with its reference
# ---------------------------------------------------------------------------

#' Text form of a table for comparison: the given columns as strings, invalid
#' UTF-8 bytes removed, missing values as "", numeric columns written the same
#' way whatever their source format ("15.0" and "15" both become "15").
#' Leading and trailing spaces are removed and "NA" counts as missing, as
#' data.table::fread() (which wrote the cubes' inputs) does by default.
chk_normalise <- function(dt, cols = CHK_COLS) {
  dt <- data.table::as.data.table(dt)
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))
  out <- data.table::data.table(.dummy = seq_len(nrow(dt)))
  for (cc in cols) {
    v <- as.character(dt[[cc]])
    v[is.na(v)] <- ""
    bad <- !stringi::stri_enc_isutf8(v)
    if (any(bad)) v[bad] <- iconv(v[bad], "UTF-8", "UTF-8", sub = "")
    Encoding(v) <- "UTF-8"
    v <- trimws(v, whitespace = "[ ]")
    v[v == "NA"] <- ""
    if (cc %in% CHK_NUM_COLS) {
      x <- suppressWarnings(as.numeric(v))
      v[!is.na(x)] <- sprintf("%.15g", x[!is.na(x)])
    }
    data.table::set(out, j = cc, value = v)
  }
  out[, .dummy := NULL]
  out[]
}

#' The reference rows as they should look after the taxonomy harmonisation
#' (scripts/harmonise_cube_taxonomy.R): rows with old-backbone (numeric)
#' species keys whose (kingdom, species) the name table resolves get the
#' table's names and keys; rows that then share cell, year and species key
#' are one row, with occurrences summed and the smallest uncertainties.
#' Input and output are normalised text tables with CHK_HARMONISED_COLS.
chk_apply_name_table <- function(ref, name_table) {
  nt <- data.table::as.data.table(name_table)[status %in% c("in update data", "accepted", "synonym")]
  nt <- unique(nt, by = c("old_kingdom", "old_species"))
  nt <- chk_normalise(nt, c("old_kingdom", "old_species", CHK_TAXON_COLS))
  r <- data.table::copy(ref)
  r[, .i := .I]
  old <- r[grepl("^[0-9]+$", specieskey), .(.i, old_kingdom = kingdom, old_species = species)]
  m <- merge(old, nt, by = c("old_kingdom", "old_species"))
  for (cc in CHK_TAXON_COLS) data.table::set(r, i = m$.i, j = cc, value = m[[cc]])
  r[, .i := NULL]
  num_min <- function(v) { x <- suppressWarnings(as.numeric(v)); if (all(is.na(x))) "" else sprintf("%.15g", min(x, na.rm = TRUE)) }
  key <- c("mgrscellcode", "year", "specieskey")
  r[, .(kingdom = kingdom[1], kingdomkey = kingdomkey[1], phylum = phylum[1], phylumkey = phylumkey[1],
        class = class[1], classkey = classkey[1], order = order[1], orderkey = orderkey[1],
        family = family[1], familykey = familykey[1], genus = genus[1], genuskey = genuskey[1],
        species = species[1],
        occurrences = sprintf("%.15g", sum(as.numeric(occurrences), na.rm = TRUE)),
        mintemporaluncertainty = num_min(mintemporaluncertainty),
        mincoordinateuncertaintyinmeters = num_min(mincoordinateuncertaintyinmeters)), by = key][, CHK_HARMONISED_COLS, with = FALSE]
}

#' Taxon group of each row: the first of the screen's subsets
#' (scripts/lib/sufficiency_criteria.R) it belongs to, else "other".
chk_group_of <- function(dt, subsets) {
  g <- rep("other", nrow(dt))
  for (s in rev(subsets)) g[subset_mask(dt, s)] <- s
  g
}

#' Compare one site cube with its reference rows.
#' `cube`, `ref`: normalised text tables with `cols`. `outline`: the site's
#' outline (sf), used to say whether extra cells lie inside the site.
#' Rows are matched on all `cols` except the higher-taxon count columns
#' (CHK_COUNT_COLS); count_cols_differ is the number of matched rows whose
#' count columns differ.
#' Returns counts and up to `n_examples` differing rows of each kind.
chk_compare_site <- function(cube, ref, cols, outline, from_year, n_examples = 20) {
  key_match <- function(cc) {
    kc <- do.call(paste, c(as.list(cube[, cc, with = FALSE]), sep = CHK_SEP))
    kr <- do.call(paste, c(as.list(ref[, cc, with = FALSE]), sep = CHK_SEP))
    tc <- data.table::data.table(k = kc)[, .(n_cube = .N, i_cube = .I[1]), by = k]
    tr <- data.table::data.table(k = kr)[, .(n_ref = .N, i_ref = .I[1]), by = k]
    m <- merge(tc, tr, by = "k", all = TRUE)
    m[is.na(n_cube), n_cube := 0L]; m[is.na(n_ref), n_ref := 0L]
    m
  }
  core <- setdiff(cols, CHK_COUNT_COLS)
  m <- key_match(core)
  same_all <- if (length(core) < length(cols)) { ma <- key_match(cols); sum(pmin(ma$n_cube, ma$n_ref)) } else sum(pmin(m$n_cube, m$n_ref))
  m[, `:=`(extra = pmax(n_cube - n_ref, 0L), missing = pmax(n_ref - n_cube, 0L))]
  ex <- m[extra > 0]; mi <- m[missing > 0]
  ex_rows <- cube[ex$i_cube]; mi_rows <- ref[mi$i_ref]
  inside <- logical(nrow(ex_rows))
  if (nrow(ex_rows)) {
    a <- chk_cells_to_sites(ex_rows$mgrscellcode, outline)$pairs
    inside <- ex_rows$mgrscellcode %in% a$mgrscellcode
  }
  y_mi <- suppressWarnings(as.numeric(mi_rows$year))
  y_ex <- suppressWarnings(as.numeric(ex_rows$year))
  list(counts = data.table::data.table(
         rows_cube = nrow(cube), rows_ref = nrow(ref), rows_same = sum(pmin(m$n_cube, m$n_ref)),
         count_cols_differ = sum(pmin(m$n_cube, m$n_ref)) - same_all,
         extra_outside_site = sum(ex$extra[!inside]), extra_inside_site = sum(ex$extra[inside]),
         missing_before = sum(mi$missing[!is.na(y_mi) & y_mi < from_year]),
         missing_from = sum(mi$missing[is.na(y_mi) | y_mi >= from_year]),
         extra_before = sum(ex$extra[!is.na(y_ex) & y_ex < from_year]),
         extra_from = sum(ex$extra[is.na(y_ex) | y_ex >= from_year])),
       examples = data.table::rbindlist(list(
         if (nrow(ex_rows)) cbind(difference = ifelse(inside, "in cube only (cell inside site)", "in cube only (cell outside site)"), ex_rows)[seq_len(min(.N, n_examples))],
         if (nrow(mi_rows)) cbind(difference = "in sources only", mi_rows)[seq_len(min(.N, n_examples))]), fill = TRUE))
}

#' Checks on a cube by itself: duplicate cell x year x species key rows,
#' missing years and cells, invalid cell codes, and which key system the rows
#' carry before and from `from_year`.
chk_cube_internal <- function(cube, from_year) {
  y <- suppressWarnings(as.numeric(cube$year))
  num_key <- grepl("^[0-9]+$", cube$specieskey)
  data.table::data.table(
    dup_cell_year_key = sum(duplicated(cube[, .(mgrscellcode, year, specieskey)])),
    missing_year = sum(is.na(y)), missing_cell = sum(cube$mgrscellcode == ""),
    invalid_cell = sum(cube$mgrscellcode != "" & !grepl(CHK_MGRS_100M, cube$mgrscellcode)),
    missing_species = sum(cube$species == ""),
    year_min = suppressWarnings(min(y, na.rm = TRUE)), year_max = suppressWarnings(max(y, na.rm = TRUE)),
    old_keys_before = sum(num_key & !is.na(y) & y < from_year), new_keys_before = sum(!num_key & !is.na(y) & y < from_year),
    old_keys_from = sum(num_key & !is.na(y) & y >= from_year), new_keys_from = sum(!num_key & !is.na(y) & y >= from_year))
}

#' Rows and occurrences per year and taxon group, cube and reference side by side.
chk_year_group_table <- function(cube, ref, subsets) {
  f <- function(d, side) {
    if (!nrow(d)) return(NULL)
    d[, .(group = chk_group_of(.SD, subsets), year = year, occ = suppressWarnings(as.numeric(occurrences))),
      .SDcols = c("kingdom", "phylum", "class", "order")][, .(rows = .N, occurrences = sum(occ, na.rm = TRUE)), by = .(year, group)][, side := side]
  }
  t <- data.table::rbindlist(list(f(cube, "cube"), f(ref, "ref")))
  if (!nrow(t)) return(data.table::data.table(year = character(), group = character(), rows_cube = integer(), rows_ref = integer(),
                                              occ_cube = numeric(), occ_ref = numeric()))
  w <- data.table::dcast(t, year + group ~ side, value.var = c("rows", "occurrences"), fill = 0)
  for (cc in c("rows_cube", "rows_ref", "occurrences_cube", "occurrences_ref")) if (!cc %in% names(w)) data.table::set(w, j = cc, value = 0)
  data.table::setnames(w, c("occurrences_cube", "occurrences_ref"), c("occ_cube", "occ_ref"))
  w[, .(year, group, rows_cube, rows_ref, occ_cube, occ_ref)]
}

# ---------------------------------------------------------------------------
# Reference rows from one byte range of a source file
# ---------------------------------------------------------------------------

chk_bucket_of <- function(site_id, n_buckets) as.integer(sub("^site_", "", site_id)) %% n_buckets

#' Read one byte range of a source file and write, for every site, the rows
#' whose cell lies in the site, to <part_dir>/bucket_NNN.tsv (site_id + the
#' 25 cube columns, as read). `role` is "before" (keep years < from_year, the
#' continental cubes) or "from" (keep years from_year to to_year, the update
#' files). `groups`: taxon subsets (lib/sufficiency_criteria.R) to keep, or
#' "all"; the classes of the rows left out are counted, so the report shows
#' what a taxon filter removed.
#' Writes <part_dir>/done.csv with the counts when the range is finished; a
#' folder without it is from an interrupted run and is renamed, not reused.
chk_reference_part <- function(path, role, part, start, end, outlines, from_year, part_dir, n_buckets = 50,
                               chunk_lines = 1e6, to_year = Inf, groups = "all") {
  done <- file.path(part_dir, "done.csv")
  if (file.exists(done)) {
    d <- data.table::fread(done)
    if (d$n_buckets != n_buckets || d$from_year != from_year || d$to_year != to_year || d$groups != paste(groups, collapse = "+"))
      stop(part_dir, " was made with other settings (buckets ", d$n_buckets, ", years ", d$from_year, "-", d$to_year,
           ", taxa ", d$groups, "): use another --out folder")
    return(d)
  }
  if (dir.exists(part_dir)) file.rename(part_dir, paste0(part_dir, "_incomplete_", format(Sys.time(), "%Y%m%d%H%M%S")))
  dir.create(part_dir, recursive = TRUE)
  st <- list(n_empty = 0, n_wrong_fields = 0, n_bad_utf8 = 0, n_invalid_code = 0, n_no_coords = 0,
             rows_in_sites = 0, rows_kept = 0, rows_other_period = 0, rows_other_taxa = 0, examples = character())
  other_classes <- list()
  t0 <- Sys.time()
  r <- chk_read_range(path, start, end, chunk_lines = chunk_lines, on_chunk = function(lines) {
    p <- chk_parse_lines(lines)
    st$n_empty <<- st$n_empty + p$n_empty; st$n_wrong_fields <<- st$n_wrong_fields + p$n_wrong_fields
    st$n_bad_utf8 <<- st$n_bad_utf8 + p$n_bad_utf8
    if (length(st$examples) < 3) st$examples <<- utils::head(c(st$examples, p$bad_examples), 3)
    if (!nrow(p$rows)) return(invisible())
    a <- chk_cells_to_sites(p$rows$mgrscellcode, outlines)
    st$n_invalid_code <<- st$n_invalid_code + a$n_invalid_code; st$n_no_coords <<- st$n_no_coords + a$n_no_coords
    if (!nrow(a$pairs)) return(invisible())
    m <- merge(a$pairs, p$rows, by = "mgrscellcode", allow.cartesian = TRUE)
    y <- suppressWarnings(as.numeric(m$year))
    keep <- if (role == "before") !is.na(y) & y < from_year else !is.na(y) & y >= from_year & y <= to_year
    in_taxa <- if ("all" %in% groups) rep(TRUE, nrow(m)) else in_subsets(m, groups)
    st$rows_in_sites <<- st$rows_in_sites + nrow(m)
    st$rows_other_period <<- st$rows_other_period + sum(!keep)
    st$rows_other_taxa <<- st$rows_other_taxa + sum(keep & !in_taxa)
    if (any(keep & !in_taxa)) other_classes[[length(other_classes) + 1]] <<- m[keep & !in_taxa, .N, by = .(class)]
    keep <- keep & in_taxa
    st$rows_kept <<- st$rows_kept + sum(keep)
    m <- m[keep, c("site_id", CHK_COLS), with = FALSE]
    if (!nrow(m)) return(invisible())
    b <- chk_bucket_of(m$site_id, n_buckets)
    for (bb in unique(b)) {
      f <- file.path(part_dir, sprintf("bucket_%03d.tsv", bb))
      data.table::fwrite(m[b == bb], f, sep = "\t", quote = FALSE, append = file.exists(f))
    }
  })
  oc <- data.table::rbindlist(other_classes)
  oc <- if (nrow(oc)) oc[, .(N = sum(N)), by = class][order(-N)] else oc
  out <- data.table::data.table(source = path, role = role, part = part, start = start, end = end,
                                n_buckets = n_buckets, from_year = from_year, to_year = to_year,
                                groups = paste(groups, collapse = "+"),
                                first = r$first, after = r$after, lines = r$lines,
                                n_empty = st$n_empty, n_wrong_fields = st$n_wrong_fields, n_bad_utf8 = st$n_bad_utf8,
                                n_invalid_code = st$n_invalid_code, n_no_coords = st$n_no_coords,
                                rows_in_sites = st$rows_in_sites, rows_kept = st$rows_kept,
                                rows_other_period = st$rows_other_period, rows_other_taxa = st$rows_other_taxa,
                                other_taxa_classes = paste(sprintf("%s:%d", oc$class, oc$N), collapse = "|"),
                                wrong_field_examples = paste(substr(st$examples, 1, 300), collapse = " || "),
                                minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  tmp <- paste0(done, ".tmp"); data.table::fwrite(out, tmp); file.rename(tmp, done)
  out
}

#' Reference rows of one bucket from all finished ranges, as a character table
#' with site_id and the cube columns.
chk_read_reference_bucket <- function(ref_dir, bucket) {
  dirs <- list.dirs(ref_dir, recursive = FALSE)
  dirs <- dirs[file.exists(file.path(dirs, "done.csv"))]
  fs <- file.path(dirs, sprintf("bucket_%03d.tsv", bucket))
  fs <- fs[file.exists(fs)]
  if (!length(fs)) return(data.table::data.table(matrix(character(), 0, length(CHK_COLS) + 1, dimnames = list(NULL, c("site_id", CHK_COLS)))))
  data.table::rbindlist(lapply(fs, data.table::fread, sep = "\t", quote = "", colClasses = "character",
                               showProgress = FALSE))
}

#' Read a site cube (all columns as text). Older cubes have extra site columns
#' in front; only the cube columns are kept.
chk_read_cube <- function(files) {
  data.table::rbindlist(lapply(files, function(f) {
    miss <- setdiff(CHK_COLS, names(data.table::fread(f, nrows = 0)))
    if (length(miss)) stop("Cube ", f, " lacks columns: ", paste(miss, collapse = ", "))
    data.table::fread(f, select = CHK_COLS, colClasses = "character", showProgress = FALSE)[, CHK_COLS, with = FALSE]
  }))
}

#' Add up "class:count|class:count" strings (classes left out by a taxon filter).
chk_sum_classes <- function(x) {
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return("")
  p <- strsplit(unlist(strsplit(x, "|", fixed = TRUE)), ":", fixed = TRUE)
  d <- data.table::data.table(class = vapply(p, function(z) paste(z[-length(z)], collapse = ":"), ""),
                              n = as.numeric(vapply(p, function(z) z[length(z)], "")))[, .(n = sum(n)), by = class][order(-n)]
  paste(sprintf("%s:%s", d$class, format(d$n, scientific = FALSE, trim = TRUE)), collapse = "|")
}
