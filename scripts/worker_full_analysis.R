#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  cat("Usage: worker_full_analysis.R <continent> <country> <site_id> <subset> <run_ts> <run_map>\n")
  quit(status = 1)
}

cont <- args[1]
country <- args[2]
site_id <- args[3]
subset <- args[4]
run_ts <- as.logical(args[5])
run_map <- as.logical(args[6])

library(b3gbi)
library(dplyr)
library(sf)

# Fast version of compute_tax_distinct_formula to speed up O(N^2) pairwise comparisons
fast_compute_tax_distinct_formula <- function(x, y) {
  temp <- names(y) %in% x$scientificName
  tax_hier_temp <- y[c(temp)]
  # Keep only valid data frames (discard NA for not-found taxa)
  is_df <- sapply(tax_hier_temp, is.data.frame)
  tax_hier_temp <- tax_hier_temp[is_df]
  
  n_spec <- length(tax_hier_temp)
  if (n_spec < 3) {
    return(NA)
  }
  L <- max(sapply(tax_hier_temp, nrow))
  all_names <- unlist(lapply(tax_hier_temp, function(df) unique(df$name)))
  counts <- table(all_names)
  shared_sum <- sum(choose(counts[counts > 1], 2))
  num_pairs <- n_spec * (n_spec - 1) / 2
  sum_of_distances <- L * num_pairs - shared_sum
  denominator <- L * num_pairs
  return(sum_of_distances / denominator)
}
assignInNamespace("compute_tax_distinct_formula", fast_compute_tax_distinct_formula, ns = "b3gbi")

# Global cache for GBIF classifications to prevent redundant queries during bootstrapping
global_classification_cache <- list()

fast_my_classification <- function(x, ...) {
  cached_names <- names(global_classification_cache)
  if (all(x %in% cached_names)) {
    return(global_classification_cache[x])
  } else {
    missing_names <- setdiff(x, cached_names)
    if (length(missing_names) > 0) {
      res <- taxize::classification(missing_names, ...)
      global_classification_cache[missing_names] <<- res
    }
    return(global_classification_cache[x])
  }
}
assignInNamespace("my_classification", fast_my_classification, ns = "b3gbi")


out_dir <- file.path("output/full_analysis", site_id)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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

data_dir <- file.path("inst/extdata", paste0("ramsar_site_data_100m_", cont))
cube_dir <- file.path(data_dir, country)
cube_files <- list.files(cube_dir, pattern = "\\.csv$", full.names = TRUE)

cube_path <- NULL
for (f in cube_files) {
  site_num <- strsplit(site_id, "_")[[1]][2] 
  if (grepl(paste0("^", site_num, "(_|\\.)"), basename(f)) || grepl(paste0("^site_", site_num, "(_|\\.)"), basename(f))) {
    cube_path <- f
    break
  }
}
if (is.null(cube_path)) quit(status = 0)

data <- tryCatch(read.csv(cube_path, stringsAsFactors = FALSE), error = function(e) NULL)
if (is.null(data) || nrow(data) < 10) quit(status = 0)
data <- data[data$year >= 1990, ]

taxa_filter <- get_taxa_filter(subset)
if (!is.null(taxa_filter)) {
  data <- tryCatch(data[taxa_filter(data), ], error = function(e) data[FALSE, ])
}

if (nrow(data) < 10) quit(status = 0)

temp_cube_file <- paste0("temp_worker_full_", site_id, "_", subset, ".csv")
write.csv(data, temp_cube_file, row.names = FALSE)

cube <- tryCatch({
  b3gbi::process_cube(temp_cube_file, cols_cellCode = "mgrscellcode", cols_speciesKey = "specieskey", cols_scientificName = "species")
}, error = function(e) NULL)

if (is.null(cube)) {
  if (file.exists(temp_cube_file)) file.remove(temp_cube_file)
  quit(status = 0)
}

indicators <- as.list(b3gbi::available_indicators)

for (ind_name in names(indicators)) {
  ind_info <- indicators[[ind_name]]
  
  if (run_ts && !is.null(ind_info$ts_wrapper) && !is.na(ind_info$ts_wrapper)) {
    out_file <- file.path(out_dir, paste0(subset, "_", ind_info$indicator_class, "_ts.rds"))
    if (!file.exists(out_file)) {
      cat(paste("Running", ind_info$ts_wrapper, "...\n"))
      res <- tryCatch({
        # Try running with CIs
        do.call(ind_info$ts_wrapper, list(data = cube, ci_type = "norm", num_bootstrap = 1000, cell_size = "auto"))
      }, error = function(e) {
        cat("CI run failed for", ind_info$ts_wrapper, ":", conditionMessage(e), "- falling back to none\n")
        # Fall back to no CIs
        tryCatch({
          do.call(ind_info$ts_wrapper, list(data = cube, ci_type = "none", cell_size = "auto"))
        }, error = function(e2) {
          cat("Fallback run failed:", conditionMessage(e2), "\n")
          NULL
        })
      })
      if (!is.null(res)) saveRDS(res, out_file)
    }
  }
}

if (file.exists(temp_cube_file)) file.remove(temp_cube_file)
quit(status = 0)
