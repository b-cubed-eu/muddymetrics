#!/usr/bin/env Rscript
# =============================================================================
# worker_temporal_only.R
#
# Recomputes ONLY the quantities that are wrong, plus a built-in regression
# test. Does NOT recompute chao2, sac_slope, precision or site area — those are
# correct in the published results and are left completely untouched.
#
# Emits per site x subset:
#   temporal_distance_cumsum  - the OLD statistic, reproduced deliberately.
#                               This is the regression test: it must match the
#                               published `temporal_distance` column. If it
#                               does, this worker is reading and processing the
#                               cubes identically to the original run, and the
#                               corrected value below can be trusted.
#   temporal_distance_annual  - the CORRECTED statistic (no cumsum).
#   total_occurrences         - sum of the `occurrences` column, i.e. exactly
#                               the numerator the original density used. Lets
#                               density be recomputed as
#                               total_occurrences / site_area_km2, using the
#                               already-correct site_area_km2 from the results.
#
# Rationale: the gates are independent, so only the broken ones need redoing.
# Recomputing the others would risk perturbing verified numbers through
# incidental differences in reimplementation.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: worker_temporal_only.R <continent> <country> <site_id> <out_file>\n")
  quit(status = 1)
}

cont <- args[1]
country <- args[2]
site_id <- args[3]
out_file <- args[4]

library(b3gbi)

SUBSETS <- c("full", "aves", "mammalia", "herps", "insects", "plants")

# Unchanged from the original worker.
get_taxa_filter <- function(subset) {
  switch(subset,
    "full" = NULL,
    "aves" = function(df) df$class == "Aves",
    "mammalia" = function(df) df$class == "Mammalia",
    "herps" = function(df) df$class %in% c("Amphibia", "Reptilia"),
    "insects" = function(df) df$class == "Insecta",
    "plants" = function(df) df$kingdom == "Plantae",
    stop("Unknown subset: ", subset)
  )
}

# Unchanged from the original worker. The bug was in what was passed in.
calculate_jaccard <- function(occ_vec, rich_vec) {
  if (length(occ_vec) < 3) return(NA_real_)

  occ_norm <- (occ_vec - min(occ_vec)) / (max(occ_vec) - min(occ_vec) + 1e-10)
  rich_norm <- (rich_vec - min(rich_vec)) / (max(rich_vec) - min(rich_vec) + 1e-10)

  occ_threshold <- mean(occ_norm)
  rich_threshold <- mean(rich_norm)

  occ_set <- which(occ_norm > occ_threshold)
  rich_set <- which(rich_norm > rich_threshold)

  intersection <- length(intersect(occ_set, rich_set))
  union_set <- length(union(occ_set, rich_set))

  if (union_set == 0) return(NA_real_)
  jaccard_sim <- intersection / union_set
  return(1 - jaccard_sim)
}

# ===================== Main =====================
processed_records <- tryCatch(read.csv(out_file, stringsAsFactors = FALSE), error = function(e) data.frame())
processed_ids <- if (nrow(processed_records) > 0 && "site_id" %in% names(processed_records)) {
  paste(processed_records$site_id, processed_records$subset, sep = "_")
} else c()

data_dir <- file.path("inst/extdata", paste0("ramsar_site_data_100m_", cont))
cube_dir <- file.path(data_dir, country)
cube_files <- list.files(cube_dir, pattern = "\\.csv$", full.names = TRUE)

cube_path <- NULL
for (f in cube_files) {
  if (grepl(site_id, basename(f))) {
    cube_path <- f
    break
  }
}
if (is.null(cube_path)) quit(status = 0)

# Identical loading and filtering to the original worker.
data <- tryCatch(read.csv(cube_path, stringsAsFactors = FALSE), error = function(e) NULL)
if (is.null(data) || nrow(data) < 10) quit(status = 0)
data <- data[data$year >= 1990, ]
if (nrow(data) < 10) quit(status = 0)

for (subset in SUBSETS) {
  if (paste(site_id, subset, sep = "_") %in% processed_ids) next

  taxa_filter <- get_taxa_filter(subset)
  if (is.null(taxa_filter)) {
    subset_data <- data
  } else {
    subset_data <- tryCatch(data[taxa_filter(data), ], error = function(e) data[FALSE, ])
  }

  if (nrow(subset_data) < 10) next

  total_occ <- sum(subset_data$occurrences, na.rm = TRUE)

  cube <- tryCatch({
    b3gbi::process_cube(subset_data,
                        cols_cellCode = "mgrscellcode",
                        cols_speciesKey = "specieskey",
                        cols_scientificName = "species")
  }, error = function(e) NULL)

  td_cumsum <- NA_real_
  td_annual <- NA_real_

  if (!is.null(cube)) {
    tot_occ_ts  <- tryCatch(as.data.frame(b3gbi::total_occ_ts(cube)[["data"]]), error = function(e) NULL)
    obs_rich_ts <- tryCatch(as.data.frame(b3gbi::obs_richness_ts(cube)[["data"]]), error = function(e) NULL)

    if (!is.null(tot_occ_ts) && !is.null(obs_rich_ts) && nrow(tot_occ_ts) >= 3) {
      merged_ts <- merge(tot_occ_ts, obs_rich_ts, by = "year")
      merged_ts <- merged_ts[order(merged_ts$year), ]
      if (nrow(merged_ts) >= 3) {
        # Regression test: exactly what the original worker did.
        td_cumsum <- calculate_jaccard(cumsum(merged_ts$diversity_val.x),
                                       cumsum(merged_ts$diversity_val.y))
        # Corrected: annual series as specified.
        td_annual <- calculate_jaccard(merged_ts$diversity_val.x,
                                       merged_ts$diversity_val.y)
      }
    }
  }

  cat(sprintf("%s,\"%s\",%s,%s,%s,%s,%s\n",
              cont, country, site_id, subset,
              ifelse(is.na(td_cumsum), "NA", td_cumsum),
              ifelse(is.na(td_annual), "NA", td_annual),
              ifelse(is.na(total_occ), "NA", total_occ)),
      file = out_file, append = TRUE)
}
quit(status = 0)
