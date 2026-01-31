# ============================================================================
# MASTER WORKFLOW: Ramsar Biodiversity Indicator Pipeline
#
# This script orchestrates the entire pipeline:
# 1. Geometry Preprocessing (Shapefile to WKT)
# 2. Data Acquisition (GBIF Downloads)
# 3. Data Processing (MGRS Splitting)
# 4. Indicator Calculation and Visualization
# ============================================================================

# --- 0. Setup ---
# Source all modular functions
lapply(list.files("R", full.names = TRUE), source)

# Load required libraries
library(b3gbi)
library(sf)
library(dplyr)

# --- 1. Geometry Preprocessing ---
# This step converts the global Ramsar shapefile into individual site WKT files.
# Run once or when the shapefile updates.
# source("scripts/process_ramsar_shapefile.R")

# --- 2. Data Acquisition ---
# Note: Requires GBIF credentials set in options() or environment variables.
# These scripts are usually run interactively or as long-running jobs.
# source("scripts/download_gbif_ramsar.R")

# --- 3. Data Processing ---
# Intersects continental MGRS cubes with Ramsar polygons.
# source("scripts/split_data_cubes.R")

# --- 4. Indicator Calculation (Example for one region) ---
# This section demonstrates how to use the modular functions for batch processing.

input_base <- "inst/extdata/ramsar_site_data_100m_asia"
output_base <- "output/ramsar_metric_results_100m/asia"
wkt_base <- "inst/extdata/ramsar_sites_wkt"

message("Starting indicator calculation for Asia...")

# Calculate Observed Richness (Time Series)
# We use plot_args to disable smoothing if data is sparse to avoid warnings.
obs_richness_results <- calc_ramsar_indicator(
  indicator = "obs_richness_ts",
  inputdir = input_base,
  maindir = output_base,
  shapefiledir = wkt_base,
  continent = "Asia",
  plot_args = list(smoothed_trend = FALSE)
)

saveRDS(obs_richness_results, file.path(output_base, "asia_obs_richness_summary.rds"))

message("Workflow complete.")
