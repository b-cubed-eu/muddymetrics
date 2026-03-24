#!/usr/bin/env Rscript
# Calculate species group specific density (1990+)
# Applies threshold per taxonomic class

library(sf)
library(units)

year_threshold <- 1990
density_threshold <- 0.25

results <- data.frame()

for(cont in c('africa','antarctica','asia','europe','northamerica','oceania','southamerica')){
  cont_dir <- paste0('output/ramsar_metric_results_100m/', cont)
  data_dir <- paste0('inst/extdata/ramsar_site_data_100m_', cont)
  wkt_dir <- paste0('inst/extdata/ramsar_sites_wkt')
  
  if(!dir.exists(cont_dir)) next
  
  files <- list.files(cont_dir, pattern='total_occ_map.RData$', recursive=TRUE)
  
  for(f in files){
    site_id <- sub('_total_occ_map.RData$', '', basename(f))
    country <- dirname(f)
    country <- basename(country)
    
    csv_path <- file.path(data_dir, country, paste0(site_id, '_data.csv'))
    if(!file.exists(csv_path)) next
    
    # Get area
    wkt_file <- file.path(wkt_dir, country, paste0(site_id, '.wkt'))
    area <- NA
    
    if(file.exists(wkt_file)){
      tryCatch({
        wkt <- paste(readLines(wkt_file), collapse='')
        poly <- st_as_sfc(wkt)
        st_crs(poly) <- 4326
        area <- as.numeric(st_area(poly)) / 1e6
      }, error = function(e) NA)
    }
    
    if(is.na(area) || area <= 0) next
    
    # Load and filter data
    tryCatch({
      dat <- read.csv(csv_path)
      dat_1990 <- dat[dat$year >= year_threshold, ]
      
      # Aggregate by class
      class_occ <- aggregate(occurrences ~ class, data = dat_1990, FUN = sum)
      
      for(i in 1:nrow(class_occ)){
        class_name <- class_occ$class[i]
        occ <- class_occ$occurrences[i]
        density <- occ / area
        
        results <- rbind(results, data.frame(
          continent = cont,
          site_id = site_id,
          country = country,
          class = class_name,
          density = density,
          passes = as.numeric(density >= density_threshold),
          area_km2 = area,
          occ_count = occ,
          stringsAsFactors = FALSE
        ))
      }
      
    }, error = function(e) NA)
  }
  
  cat('Processed', cont, ':', nrow(results[results$continent==cont,]), 'site-class combinations\n')
}

# Save full results
write.csv(results, 'output/data_sufficiency/density_by_class_1990.csv', row.names=FALSE)

cat('\n===== SPECIES GROUP DENSITY ANALYSIS =====\n')
cat('Total site-class combinations:', nrow(results), '\n')
cat('Passing (>=0.25/km2):', sum(results$passes == 1), '\n')

cat('\nBy class:\n')
class_summary <- aggregate(cbind(passes, occ_count) ~ class, data = results, FUN = sum)
class_summary$total <- table(results$class)[class_summary$class]
class_summary$pass_rate <- round(class_summary$passes / class_summary$total * 100, 1)
class_summary <- class_summary[order(-class_summary$passes), ]
print(class_summary[, c('class', 'total', 'passes', 'pass_rate')])

# Example
cat('\nExample - Algeria 1052 (Aves):\n')
print(subset(results, grepl('1052', site_id) & continent == 'africa' & class == 'Aves'))

# Save sites that have at least one passing class
sites_with_pass <- aggregate(passes ~ continent + site_id + country + area_km2, data = results, FUN = max)
cat('\nSites with at least one passing class:', sum(sites_with_pass$passes == 1), '\n')

write.csv(sites_with_pass[sites_with_pass$passes == 1, ], 
          'output/data_sufficiency/sites_for_b3gbi_by_class_1990.csv', row.names=FALSE)

cat('\nSaved to output/data_sufficiency/\n')
