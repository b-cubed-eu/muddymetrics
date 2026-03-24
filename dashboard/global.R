# Global.R - Global variables and data loading for Shiny dashboard

# Load required libraries
library(shiny)
library(leaflet)
library(dplyr)
library(ggplot2)
library(viridis)
library(plotly)
library(DT)
library(bslib)

# Data file path
summary_file <- "output/global_sufficiency_summary.csv"

# Check if data file exists
if (file.exists(summary_file)) {
  global_summary <- read.csv(summary_file, stringsAsFactors = FALSE)
} else {
  global_summary <- data.frame(
    site_id = character(),
    site_name = character(),
    country = character(),
    continent = character(),
    density_km2 = numeric(),
    data_class = character(),
    stringsAsFactors = FALSE
  )
  warning("Summary data file not found. Dashboard may not display correctly.")
}
