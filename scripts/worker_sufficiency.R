#!/usr/bin/env Rscript
# Computes all data-sufficiency criteria for one site, for each taxonomic
# subset, and writes one CSV for the site. Thresholds are NOT applied here;
# see apply_sufficiency_thresholds.R.
#
# Usage (from the repository root):
#   Rscript scripts/worker_sufficiency.R <continent> <country> <site_id> <cube_path> <out_file> [grid]
#
# grid is "100m" (default, the cube's native cells) or "1km". Only the
# criteria that use grid cells as sampling units change with the grid: Chao2
# completeness and annual coverage. Density, temporal decoupling and precision
# do not depend on it.
#
# The cube path is passed explicitly by the runner, so the worker never looks
# files up by name.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 5:6) {
  cat("Usage: worker_sufficiency.R <continent> <country> <site_id> <cube_path> <out_file> [grid]\n")
  quit(status = 1)
}
cont <- args[1]; country <- args[2]; site_id <- args[3]
cube_path <- args[4]; out_file <- args[5]
grid <- if (length(args) == 6) args[6] else "100m"

suppressPackageStartupMessages({
  library(data.table)
  library(sf)
})
source("scripts/lib/sufficiency_criteria.R")

WKT_ROOT <- "inst/extdata/ramsar_sites_wkt"
t0 <- Sys.time()

cols <- c("kingdom", "class", "species", "specieskey", "mgrscellcode", "year",
          "occurrences", "mincoordinateuncertaintyinmeters")
df_all <- as.data.frame(data.table::fread(
  cube_path, select = cols, showProgress = FALSE,
  colClasses = list(character = c("kingdom", "class", "species", "specieskey", "mgrscellcode"),
                    numeric = c("year", "occurrences", "mincoordinateuncertaintyinmeters"))))
df <- recode_grid(filter_years(df_all), grid)

wkt_files <- find_wkt(WKT_ROOT, site_id)
area <- if (length(wkt_files) >= 1) {
  tryCatch(site_area_km2(wkt_files[1]), error = function(e) NA_real_)
} else NA_real_

prec_d61 <- precision_ratio_d61(df_all, area)

rows <- lapply(SUBSETS, function(subset) {
  sub <- subset_taxon(df, subset)
  base <- list(continent = cont, country = country, site_id = site_id, subset = subset, grid = grid,
               n_rows = nrow(sub), n_occurrences = sum(sub$occurrences, na.rm = TRUE),
               site_area_km2 = area, n_wkt_matches = length(wkt_files),
               precision_ratio_d61 = prec_d61)
  if (nrow(sub) < MIN_ROWS) return(c(base, status = "too_few_rows"))

  st <- incidence_stats(sub)
  cov_tab <- annual_coverage_table(sub)
  series <- annual_series(sub)
  prec <- precision_ratio(sub, area)
  c(base,
    status = "computed",
    density_km2 = density_km2(sub, area),
    S_obs = st$S, n_cells = st$T, Q1 = st$Q1, Q2 = st$Q2, U = st$U,
    chao2 = chao2_estimate(st),
    chao2_completeness = chao2_completeness(st),
    chao2_completeness_bykey = chao2_completeness(incidence_stats(sub, "specieskey")),
    n_names_multi_key = names_with_multiple_keys(sub),
    pooled_coverage = pooled_coverage(st),
    n_years = nrow(series),
    n_years_single_cell = sum(cov_tab$n_cells == 1),
    mean_cov_b3gbi = mean_annual_coverage(cov_tab, "b3gbi"),
    mean_cov_min_cells = mean_annual_coverage(cov_tab, "min_cells"),
    mean_cov_informative = mean_annual_coverage(cov_tab, "informative"),
    n_years_informative = sum(informative_years(cov_tab)),
    temporal_cha = cha_distance(series$total_occ, series$richness),
    temporal_binary = binary_jaccard_distance(series$total_occ, series$richness),
    precision_ratio = prec$ratio,
    mean_uncertainty_area_km2 = prec$mean_area_km2,
    precision_ratio_median = prec$ratio_median,
    median_uncertainty_area_km2 = prec$median_area_km2,
    frac_occ_no_uncertainty = prec$frac_occ_missing)
})

out <- data.table::rbindlist(lapply(rows, as.data.frame), fill = TRUE)
out$coverage_deficit_b3gbi <- 1 - out$mean_cov_b3gbi
out$coverage_deficit_min_cells <- 1 - out$mean_cov_min_cells
out$coverage_deficit_informative <- 1 - out$mean_cov_informative
out$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

# Write to a temporary file and rename, so an interrupted run never leaves a
# partial file that the runner would mistake for a finished site.
tmp <- paste0(out_file, ".tmp")
data.table::fwrite(out, tmp)
invisible(file.rename(tmp, out_file))
quit(status = 0)
