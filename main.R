# ============================================================================
# MASTER WORKFLOW: Ramsar Biodiversity Indicator Pipeline
#
# This script orchestrates the entire pipeline. There are two primary paths:
#
# PATH A: Continental (Most Efficient for Thousands of Sites)
# 1. scripts/download_continental_cubes.R (Resumable 100GB+ downloads)
# 2. scripts/split_data_cubes.R           (Local cutting for Ramsar sites)
#
# PATH B: Site-by-Site (Best for testing or single-site updates)
# 1. scripts/download_gbif_ramsar_by_site.R (Direct GBIF cutting)
# ============================================================================

# --- 0. Setup ---
# Source all modular functions
lapply(list.files("R", full.names = TRUE), source)

# Load required libraries
library(b3gbi)
library(sf)
library(dplyr)

# --- 1. Indicator Calculation (Example for one region) ---
# This section demonstrates how to use the modular functions for batch processing.
# Assumes data has already been "cut" using PATH A or PATH B above.

input_base <- "inst/extdata/ramsar_site_data_100m_asia"
output_base <- "output/ramsar_metric_results_100m/asia"
wkt_base <- "inst/extdata/ramsar_sites_wkt"

message("Starting indicator calculation for Asia...")

# Calculate Observed Richness (Time Series)
# Note: ci_type = "none" is used to speed up large batch runs.
obs_richness_results <- calc_ramsar_indicator(
  indicator = "obs_richness_ts",
  inputdir = input_base,
  maindir = output_base,
  shapefiledir = wkt_base,
  continent = "Asia",
  plot_args = list(smoothed_trend = FALSE) # Disable smoothing for sparse test data
)

saveRDS(obs_richness_results, file.path(output_base, "asia_obs_richness_summary.rds"))

message("Workflow complete.")