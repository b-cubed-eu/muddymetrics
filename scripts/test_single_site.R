#!/usr/bin/env Rscript
# Test: process a single site

library(b3gbi)
library(sf)
library(units)

continent <- "asia"
country_name <- "Japan"
site_file <- "site_1200_Fujimae-higata_data.csv"

input_base <- "inst/extdata"
shapefiledir <- file.path(input_base, "ramsar_sites_wkt")
inputdir <- file.path(input_base, paste0("ramsar_site_data_100m_", continent))

sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))
site_id <- paste0(country_name, "_", sitename)

cat("Processing:", site_id, "\n")

shapefilepath <- file.path(shapefiledir, country_name, paste0(sitename, ".wkt"))
cubepath <- file.path(inputdir, country_name, site_file)

cat("Files exist:", file.exists(shapefilepath), file.exists(cubepath), "\n")

# Process
cube <- b3gbi::process_cube(cubepath, separator = ",")

total_occ_map <- b3gbi::compute_indicator_workflow(
  data = cube, type = "total_occ", dim_type = "map",
  shapefile_path = shapefilepath, ne_scale = "large", region = continent,
  include_water = TRUE, shapefile_crs = 4326, ci_type = "none"
)

area <- total_occ_map$original_bbox %>% sf::st_area() %>% units::set_units("km^2")
total_occ <- sum(total_occ_map$data$diversity_val, na.rm = TRUE)
density <- as.numeric(total_occ) / as.numeric(area)

cat("Area:", area, "\n")
cat("Occurrences:", total_occ, "\n")
cat("Density:", density, "\n")
cat("Passes:", density >= 0.25, "\n")

# Save
output_dir <- "output/density_results/asia"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

density_obj <- list(
  site_id = site_id,
  sitename = sitename,
  country = country_name,
  density_km2 = density,
  total_occurrences = as.numeric(total_occ),
  area_km2 = as.numeric(area),
  threshold = 0.25,
  passes_threshold = density >= 0.25,
  timestamp = Sys.time()
)

output_file <- file.path(output_dir, paste0(site_id, "_density.RData"))
save(density_obj, file = output_file)

cat("Saved to:", output_file, "\n")
