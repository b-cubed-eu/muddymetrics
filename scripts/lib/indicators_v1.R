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
# incidence_sc, mgrs_to_1km) and the b3gbi and iNEXT packages.

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
#' richness (all by scientific name), plus occupied cells, species x cell
#' incidences and sample coverage on the requested grid.
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
    n_incidences = vapply(inc[as.character(years)], function(z) sum(z$freq), numeric(1)),
    coverage = vapply(inc[as.character(years)], function(z) incidence_sc(z$T, z$freq), numeric(1)),
    row.names = NULL)
}

#' Hill numbers (q = 0, 1, 2) per year at a fixed sample coverage.
#' Years are kept when the species x cell matrix has more than `cutoff_length`
#' elements; b3gbi's workflow passes cutoff_length = 5 for the Hill numbers. nboot = 0 gives no intervals.
#' `years`: compute only these years (each year is estimated on its own, so the
#' result for a year does not depend on which other years are computed).
#' `seed`: with nboot > 0, the random seed for each year and q is seed + 10 * year + q,
#' so intervals can be reproduced year by year.
#' Returns one row per year and q: estimate, interval, the sample size (cells)
#' needed to reach the coverage level (t), the observed sample (T), their
#' ratio and iNEXT's method label (Rarefaction / Observed / Extrapolation).
hill_indicators <- function(df, grid = c("100m", "1km"), level = 0.95, q = c(0, 1, 2),
                            nboot = 0, conf = 0.95, cutoff_length = 5, years = NULL, seed = 20260922) {
  v <- indicator_rows(df, grid)
  inc <- yearly_incidence_freq(v)
  if (!is.null(years)) inc <- inc[names(inc) %in% as.character(years)]
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
    if (nboot > 0) set.seed(seed + 10 * as.numeric(nm) + qq)
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

# ---------------------------------------------------------------------------
# Splitting interval runs (run_indicators.R --split-hours)
# ---------------------------------------------------------------------------

#' Jobs for an interval run split by year: each site x subset combination is
#' divided into k parts, part i holding the years whose remainder divided by k
#' is i (so that early and recent years are spread over the parts). k is chosen
#' so that a part takes about `hours`, from the seconds the combination took in
#' an earlier run (`prior`, a hill_ci file; times `growth` for added data), or
#' `default_k` for combinations without an earlier time. With `size` (the
#' screen's n_rows, n_cells and S_obs per site and subset), the time of those
#' combinations is predicted instead, from a log-linear fit of the earlier times
#' on these sizes (times 1.5, to be on the safe side). Jobs are ordered by
#' their expected time, longest first. The parts together give exactly the
#' results of one run, because every year and order has its own seed.
ci_job_plan <- function(combos, prior = NULL, hours = 0.5, default_k = 8L, max_k = 40L, growth = 1.3, size = NULL) {
  cb <- unique(data.table::as.data.table(combos)[, .(site_id, subset)])
  est <- data.table::data.table(site_id = character(), subset = character(), est_hours = numeric())
  if (!is.null(prior) && nrow(prior)) {
    p <- unique(data.table::as.data.table(prior)[, .(site_id, subset, grid, level, seconds)], by = c("site_id", "subset", "grid", "level"))
    est <- p[, .(est_hours = growth * sum(seconds) / 3600), by = .(site_id, subset)]
  }
  cb <- merge(cb, est, by = c("site_id", "subset"), all.x = TRUE, sort = FALSE)
  cb[, predicted := FALSE]
  if (!is.null(size) && all(c("n_rows", "n_cells", "S_obs") %in% names(size)) && any(is.na(cb$est_hours))) {
    sz <- unique(data.table::as.data.table(size)[n_rows > 0 & n_cells > 0 & S_obs > 0, .(site_id, subset, n_rows, n_cells, S_obs)],
                 by = c("site_id", "subset"))
    known <- merge(est, sz, by = c("site_id", "subset"))
    if (nrow(known) >= 10) {
      fit <- stats::lm(log(est_hours) ~ log(n_rows) + log(n_cells) + log(S_obs), data = known)
      un <- merge(cb[is.na(est_hours), .(site_id, subset)], sz, by = c("site_id", "subset"))
      if (nrow(un)) {
        un[, pred := 1.5 * exp(stats::predict(fit, newdata = un))]
        i <- match(paste(un$site_id, un$subset), paste(cb$site_id, cb$subset))
        cb[i, `:=`(est_hours = un$pred, predicted = TRUE)]
      }
    }
  }
  cb[, k := as.integer(ifelse(is.na(est_hours), default_k, pmin(max_k, pmax(1, ceiling(est_hours / hours)))))]
  jobs <- cb[, .(part = seq_len(k) - 1L), by = .(site_id, subset, k, est_hours, predicted)]
  jobs[, years := sprintf("m%dr%d", k, part)]
  jobs[, est_part_hours := ifelse(is.na(est_hours), hours, est_hours / k)]
  jobs[order(-est_part_hours, site_id, subset, part)]
}

#' Merge the parts of a split interval run for one site into the table one
#' unsplit run would have written: rows ordered by subset (in `subset_order`),
#' grid, level, year and q; `seconds` is the sum over the parts.
merge_ci_parts <- function(files, subset_order) {
  x <- data.table::rbindlist(lapply(seq_along(files), function(i) {
    d <- data.table::fread(files[i])
    if (nrow(d) && "subset" %in% names(d)) d[, .part := i] else NULL
  }), fill = TRUE)
  if (!nrow(x)) return(data.table::data.table(site_id = character()))
  tot <- unique(x[, .(subset, grid, level, .part, seconds)], by = c("subset", "grid", "level", ".part"))[
    , .(total = sum(seconds)), by = .(subset, grid, level)]
  x <- merge(x, tot, by = c("subset", "grid", "level"), sort = FALSE)
  x[, `:=`(seconds = total, total = NULL, .part = NULL)]
  data.table::setcolorder(x, c("site_id", "subset", "grid", "level", "nboot", "seconds"))
  x[order(match(subset, subset_order), grid, level, year, q)]
}
