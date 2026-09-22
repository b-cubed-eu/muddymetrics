# Indicators for the v1 analysis (Part 2), computed from incidence frequencies.
#
# The Hill numbers use the same estimator as b3gbi::hill0_ts() / hill1_ts() /
# hill2_ts(): iNEXT::estimateD() on yearly species x cell incidence at a fixed
# sample coverage, called through b3gbi's own wrapper (b3gbi:::my_estimateD).
# The difference is that the input is the incidence *frequencies* (number of
# cells per species) instead of the dense species x cell matrix, which gives
# identical estimates without building a matrix that can exceed memory at
# large sites.
#
# Requires scripts/lib/sufficiency_criteria.R (valid_incidence_rows,
# incidence_sc) and the b3gbi and iNEXT packages.

#' MGRS 100 m cell code -> the MGRS 1 km cell containing it
#' (e.g. "32VNM641645" -> "32VNM6464"). Stops on codes that are not 100 m.
mgrs_to_1km <- function(code) {
  tail <- regmatches(code, regexpr("[0-9]{6}$", code))
  if (length(tail) != length(code)) stop("Not all cell codes are MGRS 100 m codes (6 trailing digits)")
  prefix <- substr(code, 1, nchar(code) - 6)
  paste0(prefix, substr(tail, 1, 2), substr(tail, 4, 5))
}

#' Rows used by all indicators here: taxon (scientific name) and cell present,
#' with the cell recoded to the requested grid.
indicator_rows <- function(df, grid = c("100m", "1km")) {
  grid <- match.arg(grid)
  v <- valid_incidence_rows(df)
  if (grid == "1km") v$mgrscellcode <- mgrs_to_1km(v$mgrscellcode)
  v
}

#' Per-year incidence frequencies: T (occupied cells) and the number of cells
#' occupied by each species.
yearly_incidence_freq <- function(v) {
  pairs <- unique(data.frame(year = v$year, taxon = v$taxon, cell = v$mgrscellcode,
                             stringsAsFactors = FALSE))
  lapply(split(pairs, pairs$year), function(y) {
    list(T = length(unique(y$cell)), freq = as.vector(table(y$taxon)))
  })
}

#' Simple yearly indicators: total occurrences, observed richness, cumulative
#' richness (all by scientific name), plus occupied cells and sample coverage
#' on the requested grid.
simple_indicators <- function(df, grid = c("100m", "1km")) {
  v <- indicator_rows(df, grid)
  years <- sort(unique(v$year))
  inc <- yearly_incidence_freq(v)
  first_seen <- tapply(v$year, v$taxon, min)
  data.frame(
    year = years,
    total_occ = vapply(years, function(y) sum(v$occurrences[v$year == y], na.rm = TRUE), numeric(1)),
    obs_richness = vapply(years, function(y) length(unique(v$taxon[v$year == y])), numeric(1)),
    cum_richness = cumsum(vapply(years, function(y) sum(first_seen == y), numeric(1))),
    n_cells = vapply(inc[as.character(years)], function(z) z$T, numeric(1)),
    coverage = vapply(inc[as.character(years)], function(z) incidence_sc(z$T, z$freq), numeric(1)),
    row.names = NULL)
}

#' Hill numbers (q = 0, 1, 2) per year at a fixed sample coverage.
#' Years are kept when the species x cell matrix has more than `cutoff_length`
#' elements; b3gbi's workflow passes cutoff_length = 5 for the Hill numbers. nboot = 0 gives no intervals.
#' Returns one row per year and q: estimate, interval, the sample size (cells)
#' needed to reach the coverage level (t), the observed sample (T), their
#' ratio and iNEXT's method label (Rarefaction / Observed / Extrapolation).
hill_indicators <- function(df, grid = c("100m", "1km"), level = 0.95, q = c(0, 1, 2),
                            nboot = 0, conf = 0.95, cutoff_length = 5) {
  v <- indicator_rows(df, grid)
  inc <- yearly_incidence_freq(v)
  keep <- vapply(inc, function(z) as.numeric(length(z$freq)) * z$T > cutoff_length, logical(1))
  empty <- data.frame(year = numeric(), q = numeric(), qD = numeric(), qD.LCL = numeric(),
                      qD.UCL = numeric(), t = numeric(), T = numeric(), t_over_T = numeric(),
                      method = character(), SC = numeric(), ci_status = character())
  if (!any(keep)) return(empty)
  x <- lapply(inc[keep], function(z) c(z$T, z$freq))
  # One call per year and q, so that a year where iNEXT's bootstrap fails (it can
  # stop with "NA/NaN argument" on some frequency patterns) only loses its
  # interval; its point estimate is then taken from a call without bootstrap.
  # Each q separately too: the failure is specific to q = 1 in the cases seen.
  one <- function(nm, qq) {
    run <- function(nb) suppressWarnings(b3gbi:::my_estimateD(x[nm], datatype = "incidence_freq",
      base = "coverage", level = level, q = qq, conf = conf, nboot = nb))
    e <- tryCatch(run(nboot), error = function(err) NULL)
    ci <- if (nboot == 0) "none" else "ok"
    if (is.null(e)) { e <- run(0); ci <- "failed" }
    e$ci_status <- ci
    e
  }
  est <- NULL
  if (nboot == 0) {
    # Without bootstrap, all years and orders in one call (as b3gbi does); per
    # year and order only if that call fails.
    est <- tryCatch({
      e <- suppressWarnings(b3gbi:::my_estimateD(x, datatype = "incidence_freq", base = "coverage",
                                                 level = level, q = q, conf = conf, nboot = 0))
      e$ci_status <- "none"; e }, error = function(err) NULL)
  }
  if (is.null(est)) est <- do.call(rbind, lapply(names(x), function(nm) do.call(rbind, lapply(q, function(qq) one(nm, qq)))))
  T_obs <- vapply(inc[keep], function(z) z$T, numeric(1))
  out <- data.frame(year = as.numeric(est$Assemblage), q = est$Order.q, qD = est$qD,
                    qD.LCL = est$qD.LCL, qD.UCL = est$qD.UCL, t = est$t,
                    T = unname(T_obs[as.character(est$Assemblage)]), method = est$Method,
                    SC = est$SC, ci_status = est$ci_status, stringsAsFactors = FALSE)
  out$t_over_T <- out$t / out$T
  out$qD.LCL[out$ci_status != "ok"] <- NA_real_
  out$qD.UCL[out$ci_status != "ok"] <- NA_real_
  out[order(out$year, out$q), c("year", "q", "qD", "qD.LCL", "qD.UCL", "t", "T", "t_over_T", "method", "SC", "ci_status")]
}
