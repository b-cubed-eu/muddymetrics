# Data-sufficiency criteria for site x taxon occurrence cubes.
#
# Pure functions: each takes a data frame of cube rows (one row per
# species x 100 m MGRS cell x year, or x 1 km cell after recode_grid()) and
# returns numbers. No file I/O except find_wkt() and site_area_km2(). Sourced
# by the worker and by the tests, so both always use the same code.
#
# Expected columns: kingdom, class, species, specieskey, mgrscellcode, year,
# occurrences, mincoordinateuncertaintyinmeters.
#
# Species are identified by scientific name (column `species`) by default.
# The cubes mix two GBIF key systems: numeric backbone keys on older records
# and alphanumeric keys on many records from 2023 onwards, so the same species
# can carry two different `specieskey` values. Counting by key would split it
# into two species; counting by name does not. TAXON_COL can be set to
# "specieskey" to reproduce key-based results.

# Taxon subsets screened by default. "full" (all taxa together) is still
# available on request, but it is not a consistent group over time: the
# 2024-2026 update covers only the groups below.
SUBSETS <- c("aves", "mammalia", "herps", "insects", "plants", "fish")
ALL_SUBSETS <- c("full", SUBSETS)

# The cubes use the GBIF backbone, where reptiles are split into Squamata,
# Testudines, Crocodylia and Rhynchocephalia. "Reptilia" is kept in case an
# older backbone version appears.
HERP_CLASSES <- c("Amphibia", "Reptilia", "Squamata", "Testudines",
                  "Crocodylia", "Rhynchocephalia", "Sphenodontia")

# Fish: jawless, cartilaginous and bony fishes. The old GBIF backbone gives
# ray-finned fishes no class (their orders sit directly under Chordata), so
# chordates without a class are counted as fish too, except the one non-fish
# order found without a class in the cubes (Copelata, larvaceans). The new
# backbone splits ray-finned fishes into Teleostei, Holostei, Chondrostei and
# Cladistii; these are the fish classes found in the fish download
# (0014182-260928105237408) and in the old-backbone cubes.
FISH_CLASSES <- c("Actinopterygii", "Teleostei", "Holostei", "Chondrostei", "Cladistii",
                  "Sarcopterygii", "Dipneusti", "Coelacanthi", "Elasmobranchii", "Holocephali",
                  "Chondrichthyes", "Petromyzonti", "Myxini", "Cephalaspidomorphi")
NON_FISH_ORDERS <- "Copelata"

TAXON_COL <- "species"
MIN_ROWS <- 10      # subsets with fewer cube rows are not evaluated
MIN_YEARS <- 3      # minimum years for the temporal criterion
FIRST_YEAR <- 1990
VAGUE_RATIO <- 0.10   # a record is vague if its uncertainty circle exceeds this share of the site area

# An informative year (coverage rule "informative"): at least this many occupied
# cells and species x cell incidences. In subsamples of well-sampled 1 km site-years,
# years below these limits too often show annual coverage >= 0.90 when their true
# coverage is below 0.85; at the limits that happens in under 5% of cases for every
# taxonomic group (scripts/check_informative_years.R).
INFORMATIVE_MIN_CELLS <- 20
INFORMATIVE_MIN_INCIDENCES <- 50

# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

#' Keep records from FIRST_YEAR to last_year. Rows with missing year are dropped.
filter_years <- function(df, first_year = FIRST_YEAR, last_year = Inf) {
  df[!is.na(df$year) & df$year >= first_year & df$year <= last_year, , drop = FALSE]
}

#' MGRS 100 m cell code -> the MGRS 1 km cell containing it
#' (e.g. "32VNM641645" -> "32VNM6464"). Stops on codes that are not 100 m.
mgrs_to_1km <- function(code) {
  tail <- regmatches(code, regexpr("[0-9]{6}$", code))
  if (length(tail) != length(code)) stop("Not all cell codes are MGRS 100 m codes (6 trailing digits)")
  prefix <- substr(code, 1, nchar(code) - 6)
  paste0(prefix, substr(tail, 1, 2), substr(tail, 4, 5))
}

#' Recode the cells of a cube to the requested grid. "100m" leaves the cube
#' unchanged; "1km" replaces each 100 m cell by the 1 km cell containing it.
#' Rows without a cell code are left as they are.
recode_grid <- function(df, grid = c("100m", "1km")) {
  grid <- match.arg(grid)
  if (grid == "100m") return(df)
  ok <- !is.na(df$mgrscellcode) & df$mgrscellcode != ""
  df$mgrscellcode[ok] <- mgrs_to_1km(df$mgrscellcode[ok])
  df
}

#' Rows belonging to a taxonomic subset. Uses %in% so that missing class or
#' kingdom values are excluded rather than producing NA rows (except for fish,
#' where a chordate without a class is a ray-finned fish in the old backbone).
subset_taxon <- function(df, subset) {
  df[subset_mask(df, subset), , drop = FALSE]
}

#' TRUE for the rows of df that belong to any of `subsets`.
in_subsets <- function(df, subsets = SUBSETS) {
  Reduce(`|`, lapply(subsets, function(s) subset_mask(df, s)), rep(FALSE, nrow(df)))
}

subset_mask <- function(df, subset) {
  if (subset == "fish" && (is.null(df$phylum) || is.null(df$order))) stop("The fish subset needs phylum and order columns")
  switch(subset,
    full     = rep(TRUE, nrow(df)),
    aves     = df$class %in% "Aves",
    mammalia = df$class %in% "Mammalia",
    herps    = df$class %in% HERP_CLASSES,
    insects  = df$class %in% "Insecta",
    plants   = df$kingdom %in% "Plantae",
    fish     = df$phylum %in% "Chordata" & (df$class %in% FISH_CLASSES |
                 ((is.na(df$class) | df$class == "") & !df$order %in% NON_FISH_ORDERS)),
    stop("Unknown subset: ", subset)
  )
}

#' Rows usable for incidence-based estimators: taxon and cell both present.
#' Adds a `taxon` column holding the chosen identifier.
valid_incidence_rows <- function(df, taxon_col = TAXON_COL) {
  tx <- df[[taxon_col]]
  ok <- !is.na(tx) & tx != "" & !is.na(df$mgrscellcode) & df$mgrscellcode != ""
  out <- df[ok, , drop = FALSE]
  out$taxon <- tx[ok]
  out
}

#' Diagnostic: number of species names recorded under more than one key.
names_with_multiple_keys <- function(df) {
  ok <- !is.na(df$species) & df$species != "" & !is.na(df$specieskey) & df$specieskey != ""
  k <- tapply(df$specieskey[ok], df$species[ok], function(x) length(unique(x)))
  sum(k > 1)
}

# ---------------------------------------------------------------------------
# Density
# ---------------------------------------------------------------------------

#' Occurrences per km^2 of site area.
density_km2 <- function(df, area_km2) {
  if (is.na(area_km2) || area_km2 <= 0) return(NA_real_)
  sum(df$occurrences, na.rm = TRUE) / area_km2
}

# ---------------------------------------------------------------------------
# Chao2 completeness and pooled sample coverage
# ---------------------------------------------------------------------------

#' Incidence statistics with grid cells as sampling units, years pooled.
#' Returns S (observed species), T (occupied cells), Q1, Q2 (species in
#' exactly one / two cells) and U (total species-cell incidences).
incidence_stats <- function(df, taxon_col = TAXON_COL) {
  v <- valid_incidence_rows(df, taxon_col)
  pairs <- unique(v[, c("mgrscellcode", "taxon")])
  occupancy <- table(pairs$taxon)               # cells per species
  list(S  = length(occupancy),
       T  = length(unique(pairs$mgrscellcode)),
       Q1 = sum(occupancy == 1),
       Q2 = sum(occupancy == 2),
       U  = nrow(pairs))
}

#' Bias-corrected Chao2 (Chao 1987; as in vegan::specpool):
#'   S + (T-1)/T * Q1^2 / (2 Q2)        if Q2 > 0
#'   S + (T-1)/T * Q1 (Q1-1) / 2        if Q2 = 0
#' Undefined (NA) with fewer than 2 occupied cells or no species.
chao2_estimate <- function(st) {
  if (st$S == 0 || st$T < 2) return(NA_real_)
  k <- (st$T - 1) / st$T
  if (st$Q2 > 0) st$S + k * st$Q1^2 / (2 * st$Q2)
  else           st$S + k * st$Q1 * (st$Q1 - 1) / 2
}

#' Chao2 completeness: S_obs / Chao2. NA when Chao2 is undefined.
chao2_completeness <- function(st) {
  est <- chao2_estimate(st)
  if (is.na(est) || est <= 0) return(NA_real_)
  st$S / est
}

#' Incidence-based sample coverage (Chao & Jost 2012), cells pooled over years:
#'   C = 1 - (Q1/U) * [ (T-1) Q1 / ((T-1) Q1 + 2 Q2) ]
#' NA with fewer than 2 cells, where the estimator is undefined.
pooled_coverage <- function(st) {
  if (st$S == 0 || st$T < 2) return(NA_real_)
  if (st$Q1 == 0) return(1)
  a <- (st$T - 1) * st$Q1
  1 - (st$Q1 / st$U) * (a / (a + 2 * st$Q2))
}

# ---------------------------------------------------------------------------
# Annual coverage (the "coverage deficit" criterion)
# ---------------------------------------------------------------------------

#' Incidence sample coverage from incidence frequencies, exactly as
#' iNEXT::DataInfo(datatype = "incidence_raw") computes SC (Chao & Jost 2012,
#' including iNEXT's rounding to 4 decimals). n_units = number of sampling
#' units (cells); freq = number of cells occupied by each species.
incidence_sc <- function(n_units, freq) {
  freq <- freq[freq > 0]
  U <- sum(freq); Q1 <- sum(freq == 1); Q2 <- sum(freq == 2)
  if (U == 0) return(NA_real_)
  Q0 <- if (Q2 == 0) (n_units - 1) / n_units * Q1 * (Q1 - 1) / 2 else (n_units - 1) / n_units * Q1^2 / 2 / Q2
  A <- if (Q1 > 0) n_units * Q0 / (n_units * Q0 + Q1) else 1
  round(1 - Q1 / U * A, 4)
}

#' Per-year incidence sample coverage (the quantity b3gbi::completeness_ts
#' takes from iNEXT::DataInfo), with the per-year cell count, the number of
#' species x cell incidences and the size of the species x cell matrix so that
#' different inclusion rules can be applied afterwards. Computed from incidence frequencies, so the dense species x
#' cell matrix is never built (it can exceed memory at large sites).
annual_coverage_table <- function(df, taxon_col = TAXON_COL) {
  v <- valid_incidence_rows(df, taxon_col)
  pairs <- unique(data.frame(year = v$year, taxon = v$taxon, cell = v$mgrscellcode,
                             stringsAsFactors = FALSE))
  years <- sort(unique(pairs$year))
  out <- data.frame(year = years, n_cells = NA_integer_, n_species = NA_integer_,
                    n_incidences = NA_integer_, matrix_size = NA_real_, sc = NA_real_)
  by_year <- split(pairs, pairs$year)
  for (i in seq_along(years)) {
    y <- by_year[[as.character(years[i])]]
    freq <- as.vector(table(y$taxon))
    n_cells <- length(unique(y$cell))
    out$n_cells[i] <- n_cells
    out$n_species[i] <- length(freq)
    out$n_incidences[i] <- sum(freq)
    out$matrix_size[i] <- as.numeric(length(freq)) * n_cells
    out$sc[i] <- incidence_sc(n_cells, freq)
  }
  out
}

#' Mean annual coverage under an inclusion rule.
#'   "b3gbi":     years with matrix_size > 5 (b3gbi's cutoff_length rule,
#'                which counts species x cell elements, not cells)
#'   "min_cells": years with at least 2 occupied cells and matrix_size > 5
#'   "informative": years with at least INFORMATIVE_MIN_CELLS occupied cells and
#'                INFORMATIVE_MIN_INCIDENCES species x cell incidences
mean_annual_coverage <- function(tab, rule = c("b3gbi", "min_cells", "informative")) {
  rule <- match.arg(rule)
  keep <- tab$matrix_size > 5
  if (rule == "min_cells") keep <- keep & tab$n_cells >= 2
  if (rule == "informative") keep <- informative_years(tab)
  vals <- tab$sc[keep]
  if (length(vals) == 0 || all(is.na(vals))) return(NA_real_)
  mean(vals, na.rm = TRUE)
}

#' Which rows of an annual coverage table are informative years.
informative_years <- function(tab, min_cells = INFORMATIVE_MIN_CELLS,
                              min_incidences = INFORMATIVE_MIN_INCIDENCES) {
  tab$n_cells >= min_cells & tab$n_incidences >= min_incidences
}

# ---------------------------------------------------------------------------
# Temporal decoupling of richness from effort
# ---------------------------------------------------------------------------

#' Annual total occurrences and observed richness, years with records only.
annual_series <- function(df, taxon_col = TAXON_COL) {
  v <- valid_incidence_rows(df, taxon_col)
  years <- sort(unique(v$year))
  occ  <- vapply(years, function(y) sum(v$occurrences[v$year == y], na.rm = TRUE), numeric(1))
  rich <- vapply(years, function(y) length(unique(v$taxon[v$year == y])), numeric(1))
  data.frame(year = years, total_occ = occ, richness = rich)
}

min_max <- function(x) {
  r <- range(x)
  if (r[2] == r[1]) return(rep(NA_real_, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

#' Cha (2007) extended-Jaccard (Tanimoto) distance between two series after
#' min-max scaling each to [0, 1]:
#'   d = sum (P - Q)^2 / ( sum P^2 + sum Q^2 - sum P Q )
#' 0 = identical shape (tightly coupled), 1 = maximally different.
#' NA with fewer than MIN_YEARS points or a constant series.
cha_distance <- function(p, q, min_points = MIN_YEARS) {
  if (length(p) != length(q)) stop("series lengths differ")
  if (length(p) < min_points) return(NA_real_)
  P <- min_max(p); Q <- min_max(q)
  if (anyNA(P) || anyNA(Q)) return(NA_real_)
  sum((P - Q)^2) / (sum(P^2) + sum(Q^2) - sum(P * Q))
}

#' The statistic used in the D6.1 code (kept only as a regression check):
#' binary Jaccard distance between the sets of above-mean years.
binary_jaccard_distance <- function(p, q, min_points = MIN_YEARS) {
  if (length(p) < min_points) return(NA_real_)
  pn <- (p - min(p)) / (max(p) - min(p) + 1e-10)
  qn <- (q - min(q)) / (max(q) - min(q) + 1e-10)
  a <- which(pn > mean(pn)); b <- which(qn > mean(qn))
  u <- length(union(a, b))
  if (u == 0) return(NA_real_)
  1 - length(intersect(a, b)) / u
}

#' Elbow threshold across all combinations (as in D6.1): sort distances,
#' scale rank to x in [0,1] and distance to y in [0,1], and take the distance
#' at maximum perpendicular deviation from the 1:1 line.
elbow_threshold <- function(distances) {
  d <- sort(distances[!is.na(distances)])
  n <- length(d)
  if (n < 10) stop("fewer than 10 distances; elbow undefined")
  x <- seq_len(n) / n
  y <- (d - d[1]) / (d[n] - d[1])
  d[which.max(abs(x - y) / sqrt(2))]
}

# ---------------------------------------------------------------------------
# Spatial precision ratio
# ---------------------------------------------------------------------------

#' Occurrence-weighted mean uncertainty area (pi r^2, km^2) over rows that
#' report an uncertainty, divided by site area. The cubes only carry the
#' minimum uncertainty per row, so this is a lower bound on true imprecision.
#' Also returns the share of occurrences with no uncertainty value.
precision_ratio <- function(df, area_km2) {
  u <- suppressWarnings(as.numeric(df$mincoordinateuncertaintyinmeters))
  w <- df$occurrences
  has <- !is.na(u) & !is.na(w) & w > 0
  total_w <- sum(w[!is.na(w) & w > 0])
  frac_missing <- if (total_w > 0) 1 - sum(w[has]) / total_w else NA_real_
  if (!any(has) || is.na(area_km2) || area_km2 <= 0) {
    return(list(ratio = NA_real_, mean_area_km2 = NA_real_, frac_occ_missing = frac_missing,
                ratio_median = NA_real_, median_area_km2 = NA_real_))
  }
  mean_area <- sum(w[has] * pi * (u[has] / 1000)^2) / sum(w[has])
  median_area <- pi * (weighted_median(u[has], w[has]) / 1000)^2
  list(ratio = mean_area / area_km2, mean_area_km2 = mean_area, frac_occ_missing = frac_missing,
       ratio_median = median_area / area_km2, median_area_km2 = median_area)
}

#' Rows whose uncertainty circle (pi r^2, from mincoordinateuncertaintyinmeters)
#' is larger than `max_ratio` times the site area: records too vague to place
#' within the site (the same ratio as the precision criterion, applied to single
#' records). Rows without an uncertainty value are not counted as vague.
vague_rows <- function(df, area_km2, max_ratio = VAGUE_RATIO) {
  u <- suppressWarnings(as.numeric(df$mincoordinateuncertaintyinmeters))
  if (is.na(area_km2) || area_km2 <= 0) return(rep(FALSE, nrow(df)))
  !is.na(u) & pi * (u / 1000)^2 > max_ratio * area_km2
}

#' Share of the rows and of the occurrences of a subset that are vague, among
#' those with an uncertainty value.
vague_share <- function(df, area_km2, max_ratio = VAGUE_RATIO) {
  u <- suppressWarnings(as.numeric(df$mincoordinateuncertaintyinmeters))
  has <- !is.na(u); v <- vague_rows(df, area_km2, max_ratio); w <- df$occurrences
  if (!any(has) || is.na(area_km2) || area_km2 <= 0) return(list(rows = NA_real_, occ = NA_real_))
  list(rows = sum(v) / sum(has), occ = sum(w[v], na.rm = TRUE) / sum(w[has], na.rm = TRUE))
}

#' Lower weighted median: the smallest x at which the cumulative weight
#' reaches half the total weight.
weighted_median <- function(x, w) {
  o <- order(x); x <- x[o]; w <- w[o]
  x[which(cumsum(w) >= sum(w) / 2)[1]]
}

#' D6.1's definition, kept only as a regression check: unweighted mean radius
#' over all rows of the whole cube (all taxa, all years), then pi r^2.
precision_ratio_d61 <- function(df_all, area_km2) {
  u <- suppressWarnings(as.numeric(df_all$mincoordinateuncertaintyinmeters))
  r <- mean(u, na.rm = TRUE)
  if (is.nan(r) || is.na(area_km2) || area_km2 <= 0) return(NA_real_)
  pi * (r / 1000)^2 / area_km2
}

# ---------------------------------------------------------------------------
# Site geometry
# ---------------------------------------------------------------------------

#' The WKT file for a site: exact match on "site_<n>" followed by "_", "." or
#' end of name, searched recursively. Returns all matches so callers can
#' detect ambiguity.
find_wkt <- function(wkt_root, site_id) {
  list.files(wkt_root, pattern = paste0("^", site_id, "(_|\\.|$)"),
             recursive = TRUE, full.names = TRUE)
}

#' Geodesic area in km^2 (WGS84 ellipsoid via lwgeom, s2 off), after
#' st_make_valid. Same method as the D6.1 precision step.
site_area_km2 <- function(wkt_path) {
  old <- suppressMessages(sf::sf_use_s2(FALSE)); on.exit(suppressMessages(sf::sf_use_s2(old)), add = TRUE)
  wkt <- paste(readLines(wkt_path, warn = FALSE), collapse = " ")
  poly <- sf::st_make_valid(sf::st_as_sfc(wkt, crs = 4326))
  as.numeric(sum(sf::st_area(poly))) / 1e6
}

# ---------------------------------------------------------------------------
# Species key harmonisation (for anything that counts species by key, e.g. b3gbi)
# ---------------------------------------------------------------------------

#' Give every species name a single `specieskey`, so that tools which count
#' species by key (including b3gbi 1.0.0) do not count one species twice when
#' it appears under both an old numeric backbone key and a new alphanumeric key.
#' Rule: per name, use the most frequent numeric key if there is one, otherwise
#' the most frequent key. Rows without a name are left unchanged.
#' Rows that end up with the same cell, year and key are then merged (see
#' merge_duplicate_cube_rows), because b3gbi::process_cube() drops identical
#' rows and would otherwise lose occurrences.
harmonise_species_keys <- function(df) {
  ok <- !is.na(df$species) & df$species != "" & !is.na(df$specieskey) & df$specieskey != ""
  pick <- function(keys) {
    num <- keys[grepl("^[0-9]+$", keys)]
    pool <- if (length(num)) num else keys
    tab <- sort(table(pool), decreasing = TRUE)
    names(tab)[1]
  }
  chosen <- tapply(df$specieskey[ok], df$species[ok], pick)
  df$specieskey[ok] <- unname(chosen[df$species[ok]])
  merge_duplicate_cube_rows(df)
}

#' Merge rows sharing cell, year, species key and name: occurrences are summed,
#' minimum uncertainties take the minimum, other columns keep the first row's
#' value. Row order of the first occurrence of each group is kept.
merge_duplicate_cube_rows <- function(df) {
  if (!all(c("mgrscellcode", "year", "specieskey", "species", "occurrences") %in% names(df))) return(df)
  k <- paste(df$mgrscellcode, df$year, df$specieskey, df$species, sep = "\r")
  if (!anyDuplicated(k)) return(df)
  g <- match(k, unique(k))
  first <- !duplicated(g)
  out <- df[first, , drop = FALSE]
  gi <- g[first]
  out$occurrences <- as.vector(tapply(df$occurrences, g, sum))[gi]
  for (col in intersect(c("mincoordinateuncertaintyinmeters", "mintemporaluncertainty"), names(df))) {
    v <- suppressWarnings(as.numeric(df[[col]]))
    m <- tapply(v, g, function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE))
    out[[col]] <- as.vector(m)[gi]
  }
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# Site cubes
# ---------------------------------------------------------------------------

#' One cube per site_id under inst/extdata/ramsar_site_data_100m_<continent>/<country>/,
#' the largest file where a site appears more than once (the D6.1 rule).
discover_site_cubes <- function(root = "inst/extdata",
                                continents = c("africa", "antarctica", "asia", "europe",
                                               "northamerica", "oceania", "southamerica")) {
  found <- list()
  for (cont in continents) {
    data_dir <- file.path(root, paste0("ramsar_site_data_100m_", cont))
    if (!dir.exists(data_dir)) next
    for (country in list.dirs(data_dir, full.names = FALSE, recursive = FALSE)) {
      for (cf in list.files(file.path(data_dir, country), pattern = "\\.csv$", full.names = TRUE)) {
        fn <- basename(cf)
        m <- regexpr("^site_[0-9]+", fn)
        site_id <- if (m != -1) regmatches(fn, m) else sub("(_[a-zA-Z0-9_]+)?_data\\.csv$", "", fn)
        found[[length(found) + 1]] <- data.frame(cont = cont, country = country, site_id = site_id,
                                                 file_path = cf, size = file.info(cf)$size,
                                                 stringsAsFactors = FALSE)
      }
    }
  }
  all <- do.call(rbind, found)
  all$order <- seq_len(nrow(all))
  x <- all[order(all$site_id, -all$size, all$order), ]
  x <- x[!duplicated(x$site_id), ]
  x <- x[order(x$order), ]
  x$order <- NULL
  rownames(x) <- NULL
  x
}
