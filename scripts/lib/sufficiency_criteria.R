# Data-sufficiency criteria for site x taxon occurrence cubes.
#
# Pure functions: each takes a data frame of cube rows (one row per
# species x 100 m MGRS cell x year) and returns numbers. No file I/O except
# find_wkt() and site_area_km2(). Sourced by the worker and by the tests, so
# both always use the same code.
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

SUBSETS <- c("full", "aves", "mammalia", "herps", "insects", "plants")

# The cubes use the GBIF backbone, where reptiles are split into Squamata,
# Testudines, Crocodylia and Rhynchocephalia. "Reptilia" is kept in case an
# older backbone version appears.
HERP_CLASSES <- c("Amphibia", "Reptilia", "Squamata", "Testudines",
                  "Crocodylia", "Rhynchocephalia")

TAXON_COL <- "species"
MIN_ROWS <- 10      # subsets with fewer cube rows are not evaluated
MIN_YEARS <- 3      # minimum years for the temporal criterion
FIRST_YEAR <- 1990

# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

#' Keep records from FIRST_YEAR onwards. Rows with missing year are dropped.
filter_years <- function(df, first_year = FIRST_YEAR) {
  df[!is.na(df$year) & df$year >= first_year, , drop = FALSE]
}

#' Rows belonging to a taxonomic subset. Uses %in% so that missing class or
#' kingdom values are excluded rather than producing NA rows.
subset_taxon <- function(df, subset) {
  keep <- switch(subset,
    full     = rep(TRUE, nrow(df)),
    aves     = df$class %in% "Aves",
    mammalia = df$class %in% "Mammalia",
    herps    = df$class %in% HERP_CLASSES,
    insects  = df$class %in% "Insecta",
    plants   = df$kingdom %in% "Plantae",
    stop("Unknown subset: ", subset)
  )
  df[keep, , drop = FALSE]
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

#' Incidence statistics with 100 m cells as sampling units, years pooled.
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

#' Per-year incidence sample coverage via iNEXT::DataInfo (as used by
#' b3gbi::completeness_ts), with the per-year cell count and matrix size so
#' that different inclusion rules can be applied afterwards.
annual_coverage_table <- function(df, taxon_col = TAXON_COL) {
  v <- valid_incidence_rows(df, taxon_col)
  years <- sort(unique(v$year))
  out <- data.frame(year = years, n_cells = NA_integer_, n_species = NA_integer_,
                    matrix_size = NA_integer_, sc = NA_real_)
  for (i in seq_along(years)) {
    y <- v[v$year == years[i], ]
    pairs <- unique(y[, c("taxon", "mgrscellcode")])
    m <- table(pairs$taxon, pairs$mgrscellcode)
    m <- (unclass(m) > 0) * 1L
    out$n_cells[i] <- ncol(m)
    out$n_species[i] <- nrow(m)
    out$matrix_size[i] <- length(m)
    sc <- tryCatch(
      suppressWarnings(iNEXT::DataInfo(list(a = m), datatype = "incidence_raw"))$SC,
      error = function(e) NA_real_)
    out$sc[i] <- as.numeric(sc)
  }
  out
}

#' Mean annual coverage under an inclusion rule.
#'   "b3gbi":     years with matrix_size > 5 (b3gbi's cutoff_length rule,
#'                which counts species x cell elements, not cells)
#'   "min_cells": years with at least 2 occupied cells and matrix_size > 5
mean_annual_coverage <- function(tab, rule = c("b3gbi", "min_cells")) {
  rule <- match.arg(rule)
  keep <- tab$matrix_size > 5
  if (rule == "min_cells") keep <- keep & tab$n_cells >= 2
  vals <- tab$sc[keep]
  if (length(vals) == 0 || all(is.na(vals))) return(NA_real_)
  mean(vals, na.rm = TRUE)
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
    return(list(ratio = NA_real_, mean_area_km2 = NA_real_, frac_occ_missing = frac_missing))
  }
  mean_area <- sum(w[has] * pi * (u[has] / 1000)^2) / sum(w[has])
  list(ratio = mean_area / area_km2, mean_area_km2 = mean_area, frac_occ_missing = frac_missing)
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
  df
}
