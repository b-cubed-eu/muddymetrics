#!/usr/bin/env Rscript
# How much information does a year need before its sample coverage and
# coverage-standardised Hill numbers can be trusted? This check sets
# INFORMATIVE_MIN_CELLS and INFORMATIVE_MIN_INCIDENCES in
# scripts/lib/sufficiency_criteria.R.
#
# Usage (from the repository root, after run_indicators.R in point mode):
#   Rscript scripts/check_informative_years.R
#
# Output: output/v1_checks/informative_years/ (reference site-years, one row per
# subsample for coverage and for the Hill numbers, and summary tables).
#
# Well-sampled site-years (>= 100 occupied 1 km cells) are treated as the
# "population". From each, k cells are drawn at random without replacement,
# many times, and the estimates from the k-cell subsample are compared with
# what the full year implies:
#   coverage: estimated sample coverage of the subsample (Chao & Jost 2012,
#             as in iNEXT) vs its true coverage relative to the full year,
#             sum of full-year incidences of the species the subsample found /
#             sum of all full-year incidences;
#   Hill numbers q = 0, 1 at 90% coverage (b3gbi:::my_estimateD, as in the
#             indicators) from the subsample vs from the full year.
# Each subsample also records the information it holds (species, incidences,
# species in one / two / two or more cells), so a minimum-information rule can
# be read off the results.

suppressPackageStartupMessages({ library(data.table); library(b3gbi) })
source("scripts/lib/sufficiency_criteria.R")

set.seed(20261006)
K <- c(2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40, 50)
R_COV <- 50; R_HILL <- 20; LEVEL <- 0.9
N_PER_SUBSET <- 50; MIN_T <- 100
OUT <- "output/v1_checks/informative_years"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

simple <- fread("output/indicators_v1/simple.csv")[grid == "1km" & n_cells >= MIN_T & coverage >= 0.95]   # reference years must themselves be well covered
ref <- simple[, .SD[sample(.N, 1)], by = .(subset, site_id)]             # one year per site
ref <- ref[, .SD[sample(.N, min(.N, N_PER_SUBSET))], by = subset]
fwrite(ref[, .(site_id, subset, year, n_cells, obs_richness)], file.path(OUT, "reference_site_years.csv"))
cat("Reference site-years:", nrow(ref), "\n"); print(ref[, .N, by = subset])

cubes <- discover_site_cubes()
cols <- c("kingdom", "class", "species", "specieskey", "mgrscellcode", "year", "occurrences")

one_site_year <- function(d, sub, yr) {
  v <- valid_incidence_rows(recode_grid(subset_taxon(d[d$year == yr, , drop = FALSE], sub), "1km"))
  p <- unique(data.table(taxon = v$taxon, cell = v$mgrscellcode))
  cells <- split(p$taxon, p$cell)
  Tn <- length(cells)
  Yfull <- table(p$taxon)
  tot <- sum(Yfull)
  rows <- list(); hill_in <- list(full = c(Tn, as.vector(Yfull))); hill_meta <- list()
  for (k in K[K <= Tn / 2]) for (r in seq_len(R_COV)) {
    pick <- cells[sample.int(Tn, k)]
    f <- table(unlist(pick, use.names = FALSE))
    fr <- as.vector(f)
    rows[[length(rows) + 1]] <- data.table(k = k, rep = r, S = length(fr), U = sum(fr),
      Q1 = sum(fr == 1), Q2 = sum(fr == 2), S2plus = sum(fr >= 2),
      C_hat = incidence_sc(k, fr), C_true = sum(Yfull[names(f)]) / tot)
    if (r <= R_HILL && length(fr) * k > 5) {
      nm <- paste(k, r, sep = "_")
      hill_in[[nm]] <- c(k, fr); hill_meta[[nm]] <- c(k, r)
    }
  }
  cov <- rbindlist(rows)
  est <- tryCatch(suppressWarnings(b3gbi:::my_estimateD(hill_in, datatype = "incidence_freq",
                    base = "coverage", level = LEVEL, q = c(0, 1), nboot = 0)), error = function(e) NULL)
  if (is.null(est)) {   # fall back to one call per subsample
    est <- rbindlist(lapply(names(hill_in), function(nm) tryCatch(as.data.table(suppressWarnings(
      b3gbi:::my_estimateD(hill_in[nm], datatype = "incidence_freq", base = "coverage",
                           level = LEVEL, q = c(0, 1), nboot = 0))), error = function(e) NULL)), fill = TRUE)
  }
  est <- as.data.table(est)[, .(assemblage = as.character(Assemblage), q = Order.q, qD, t, SC, method = Method)]
  full <- est[assemblage == "full", .(q, qD_full = qD, SC_full = SC)]
  sub_est <- est[assemblage != "full"]
  sub_est[, c("k", "rep") := tstrsplit(assemblage, "_", type.convert = TRUE)]
  hill <- merge(sub_est, full, by = "q")
  list(cov = cov, hill = hill[, .(k, rep, q, qD, qD_full, t, SC, SC_full, method)], T_full = Tn,
       S_full = length(Yfull), C_full = incidence_sc(Tn, as.vector(Yfull)))
}

res_cov <- list(); res_hill <- list()
for (s in unique(ref$site_id)) {
  path <- cubes$file_path[cubes$site_id == s]
  d <- as.data.frame(fread(path, select = cols, showProgress = FALSE,
    colClasses = list(character = c("kingdom", "class", "species", "specieskey", "mgrscellcode"))))
  d <- filter_years(d)
  for (i in which(ref$site_id == s)) {
    sub <- ref$subset[i]; yr <- ref$year[i]
    o <- tryCatch(one_site_year(d, sub, yr), error = function(e) { message(s, " ", sub, " ", yr, ": ", conditionMessage(e)); NULL })
    if (is.null(o)) next
    tag <- data.table(site_id = s, subset = sub, year = yr, T_full = o$T_full, S_full = o$S_full, C_full = o$C_full)
    res_cov[[length(res_cov) + 1]] <- cbind(tag, o$cov)
    res_hill[[length(res_hill) + 1]] <- cbind(tag, o$hill)
  }
  cat(sprintf("%s done (%d of %d site-years)\n", s, length(res_cov), nrow(ref)))
}
fwrite(rbindlist(res_cov), file.path(OUT, "subsample_coverage.csv"))
fwrite(rbindlist(res_hill), file.path(OUT, "subsample_hill.csv"))

# Summaries. A "false pass" is a subsample whose estimated coverage is >= 0.90
# while its true coverage is below 0.85: the error that lets a thin year pass
# the coverage criterion.
cv <- rbindlist(res_cov)
by_k <- cv[C_hat >= 0.9, .(n = .N, false_pass = mean(C_true < 0.85),
                           false_pass_given_50_incidences = mean(C_true[U >= 50] < 0.85)), by = k][order(k)]
inf <- cv[C_hat >= 0.9 & k >= INFORMATIVE_MIN_CELLS & k <= INFORMATIVE_MIN_CELLS + 5 & U >= INFORMATIVE_MIN_INCIDENCES,
          .(n = .N, false_pass = mean(C_true < 0.85)), by = subset]
hl <- merge(rbindlist(res_hill), cv[, .(site_id, subset, k, rep, U)], by = c("site_id", "subset", "k", "rep"))
hl <- hl[abs(SC - LEVEL) < 1e-3]
hl[, extrapolation := cut(t / k, c(0, 1, 1.5, 2, 3, Inf), labels = c("none", "1-1.5x", "1.5-2x", "2-3x", ">3x"))]
hill_ext <- hl[, .(n = .N, median_abs_rel_error = median(abs(qD / qD_full - 1)),
                   within_20pct = mean(abs(qD / qD_full - 1) <= 0.2)), by = .(q, extrapolation)][order(q, extrapolation)]
fwrite(by_k, file.path(OUT, "false_pass_by_cells.csv"))
fwrite(inf, file.path(OUT, "false_pass_at_informative_limit_by_taxon.csv"))
fwrite(hill_ext, file.path(OUT, "hill_error_by_extrapolation.csv"))
print(by_k); print(inf); print(hill_ext)
cat("Finished.\n")
