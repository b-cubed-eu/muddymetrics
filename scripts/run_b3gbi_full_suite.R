#!/usr/bin/env Rscript
# Calculate full B3GBI suite for sites passing data sufficiency threshold

library(b3gbi)
library(sf)
library(ggplot2)
library(tools)

# Read passing sites
passing_sites <- read.csv('output/data_sufficiency/sites_for_b3gbi_wkt.csv')

cat("Total passing sites:", nrow(passing_sites), "\n")

# Define all indicators in the B3GBI suite
indicators <- c(
  "total_occ",
  "obs_richness", 
  "cum_richness",
  "hill0",
  "hill1",
  "hill2",
  "ab_rarity",
  "area_rarity",
  "newness",
  "occ_density",
  "pielou_evenness",
  "williams_evenness",
  "tax_distinct"
)

# Base directories
shapefiledir <- "inst/extdata/ramsar_sites_wkt"
output_base <- "output/ramsar_metric_results_100m"

# Process by continent
processed <- 0
errors <- 0

for(cont in unique(passing_sites$continent)){
  cont_sites <- passing_sites[passing_sites$continent == cont,]
  inputdir <- paste0("inst/extdata/ramsar_site_data_100m_", cont)
  
  if(!dir.exists(inputdir)){
    cat("Skipping", cont, "- no input directory\n")
    next
  }
  
  cat("\nProcessing continent:", cont, "-", nrow(cont_sites), "sites\n")
  
  # Get unique countries in this continent
  countries <- unique(cont_sites$country)
  
  for(country in countries){
    country_sites <- cont_sites[cont_sites$country == country,]
    country_input_dir <- file.path(inputdir, country)
    
    if(!dir.exists(country_input_dir)){
      cat("  Skipping", country, "- no input directory\n")
      next
    }
    
    cat("  Processing", country, "-", nrow(country_sites), "sites\n")
    
    for(idx in 1:nrow(country_sites)){
      site_row <- country_sites[idx,]
      sitename <- site_row$site_id
      
      # Expected file paths
      cubepath <- file.path(country_input_dir, paste0(sitename, "_data.csv"))
      shapefilepath <- file.path(shapefiledir, country, paste0(sitename, ".wkt"))
      
      if(!file.exists(cubepath)){
        cat("    Missing cube:", sitename, "\n")
        errors <- errors + 1
        next
      }
      if(!file.exists(shapefilepath)){
        cat("    Missing WKT:", sitename, "\n")
        errors <- errors + 1
        next
      }
      
      # Create output directory
      output_dir <- file.path(output_base, cont, country)
      if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
      
      # Process each indicator
      for(ind in indicators){
        for(dim_type in c("map", "ts")){
          indicator_name <- paste0(ind, "_", dim_type)
          output_file <- file.path(output_dir, paste0(sitename, "_", indicator_name, ".RData"))
          
          # Skip if already exists
          if(file.exists(output_file)){
            next
          }
          
          tryCatch({
            # Process cube
            cube <- b3gbi::process_cube(cubepath, separator = ",")
            
            # Calculate indicator
            result <- b3gbi::compute_indicator_workflow(
              data = cube,
              type = ind,
              dim_type = dim_type,
              shapefile_path = shapefilepath,
              ne_scale = "large",
              region = cont,
              include_water = TRUE,
              shapefile_crs = 4326,
              ci_type = "none"
            )
            
            # Save result
            saveRDS(result, file = output_file)
            
            # Generate plot
            if(dim_type == "ts"){
              plot <- b3gbi::plot_ts(result, x_expand = 0.1, y_expand = 0.1)
            } else {
              plot <- b3gbi::plot_map(result, layers = "lakes")
            }
            
            ggplot2::ggsave(
              filename = file.path(output_dir, paste0(sitename, "_", indicator_name, ".png")),
              plot = plot,
              device = "png",
              height = 4000,
              width = 4000,
              units = "px"
            )
            
          }, error = function(e){
            cat("      Error for", indicator_name, ":", conditionMessage(e), "\n")
          })
        }
      }
      
      processed <- processed + 1
      if(processed %% 10 == 0) cat("    Processed:", processed, "sites\n")
    }
  }
}

cat("\n===== COMPLETE =====\n")
cat("Total sites processed:", processed, "\n")
cat("Errors:", errors, "\n")
