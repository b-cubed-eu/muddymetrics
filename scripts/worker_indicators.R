#!/usr/bin/env Rscript
# Computes the v1 Part 2 indicators for one site, for every taxonomic subset
# the sufficiency screen evaluated.
#
# Usage (from the repository root):
#   Rscript scripts/worker_indicators.R <site_id> <cube_path> <out_prefix> <subsets> <mode> <nboot> [grids] [years]
#     subsets  comma-separated subsets (e.g. "aves,herps,fish")
#     mode     "point": all indicators without intervals, for every subset listed
#              "ci":    Hill numbers with bootstrap intervals only, for the subsets listed
#     nboot    bootstrap replicates (mode "ci")
#     grids    optional, comma-separated grids (default "100m,1km")
#     years    optional, mode "ci" only: first and last year ("2024-2025"); only these
#              years get Hill numbers and intervals, to add years to an earlier run;
#              or "m<K>r<i>" (e.g. "m8r3"): the years whose remainder divided by K is i,
#              for splitting one combination over several processes (run_indicators.R
#              --split-hours). Each year and order has its own random seed, so the parts
#              together give exactly the results of one run over all years.
#
# Indicators (scripts/lib/indicators_v1.R): total occurrences, observed and
# cumulative richness, occupied cells and sample coverage, and Hill numbers
# q = 0, 1, 2 at 90% and 95% sample coverage, each on the 100 m and 1 km MGRS
# grids. Species are identified by scientific name. MAGV is computed from the
# annual total occurrences.
#
# Mode "point" writes <out_prefix>_simple.csv, _hill.csv and _magv.csv; mode "ci"
# writes <out_prefix>_hill_ci_<grids>_b<nboot>.csv (with _y<years> when years
# are given). Files are written under temporary names and
# renamed at the end, so a finished site always has complete files.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 6:8) {
  cat("Usage: worker_indicators.R <site_id> <cube_path> <out_prefix> <subsets> <mode> <nboot> [grids] [years]\n")
  quit(status = 1)
}
site_id <- args[1]; cube_path <- args[2]; out_prefix <- args[3]
subsets <- strsplit(args[4], ",")[[1]]
mode <- match.arg(args[5], c("point", "ci"))
nboot <- as.integer(args[6])

suppressPackageStartupMessages({
  library(data.table)
  library(b3gbi)
})
source("scripts/lib/sufficiency_criteria.R")
source("scripts/lib/indicators_v1.R")
source("scripts/lib/magv.R")

GRIDS <- if (length(args) >= 7) strsplit(args[7], ",")[[1]] else c("100m", "1km")
YEARS_ARG <- if (length(args) == 8) args[8] else ""
YEARS <- NULL
if (nzchar(YEARS_ARG)) {
  if (mode != "ci") stop("years can only be given in mode \"ci\" (cumulative richness and MAGV need all years)")
  if (grepl("^m[0-9]+r[0-9]+$", YEARS_ARG)) {
    YEARS_MOD <- as.numeric(strsplit(sub("^m", "", YEARS_ARG), "r")[[1]])   # K, i: chosen below from the cube's years
  } else {
    yy <- as.numeric(strsplit(YEARS_ARG, "-")[[1]])
    YEARS <- seq(yy[1], yy[length(yy)])
  }
}
stopifnot(all(GRIDS %in% c("100m", "1km")))
LEVELS <- c(0.90, 0.95)

cols <- c("kingdom", "phylum", "class", "order", "species", "specieskey", "mgrscellcode", "year", "occurrences")
df_all <- as.data.frame(fread(cube_path, select = cols, showProgress = FALSE,
  colClasses = list(character = c("kingdom", "phylum", "class", "order", "species", "specieskey", "mgrscellcode"),
                    numeric = c("year", "occurrences"))))
df_all <- filter_years(df_all)
if (exists("YEARS_MOD")) {
  yr <- sort(unique(df_all$year))
  YEARS <- yr[yr %% YEARS_MOD[1] == YEARS_MOD[2]]
}

simple <- hill <- magv <- list()
for (sub in subsets) {
  d <- subset_taxon(df_all, sub)
  if (nrow(d) < MIN_ROWS) next
  for (g in GRIDS) {
    if (mode == "point") {
      s <- simple_indicators(d, g)
      simple[[length(simple) + 1]] <- cbind(data.frame(site_id = site_id, subset = sub, grid = g), s)
      if (g == "100m") {
        m <- magv_series(s$year, s$total_occ)
        magv[[length(magv) + 1]] <- cbind(data.frame(site_id = site_id, subset = sub), m)
      }
    }
    for (lv in LEVELS) {
      nb <- if (mode == "ci") nboot else 0L
      t0 <- Sys.time()
      h <- hill_indicators(d, g, level = lv, nboot = nb, years = YEARS)
      if (nrow(h)) hill[[length(hill) + 1]] <- cbind(data.frame(site_id = site_id, subset = sub, grid = g,
        level = lv, nboot = nb, seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))), h)
    }
  }
}

write_tmp <- function(x, suffix) {
  tmp <- paste0(out_prefix, suffix, ".tmp")
  fwrite(if (length(x)) rbindlist(x, fill = TRUE) else data.table(site_id = character()), tmp)
  tmp
}
tmps <- if (mode == "point") {
  c(write_tmp(simple, "_simple.csv"), write_tmp(hill, "_hill.csv"), write_tmp(magv, "_magv.csv"))
} else write_tmp(hill, paste0("_hill_ci_", paste(GRIDS, collapse = "_"), "_b", nboot,
                              if (nzchar(YEARS_ARG)) paste0("_y", YEARS_ARG) else "", ".csv"))
invisible(file.rename(tmps, sub("\\.tmp$", "", tmps)))
quit(status = 0)
