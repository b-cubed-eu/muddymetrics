#!/usr/bin/env Rscript
library(sf)
library(units)

year_threshold <- 1990

results <- data.frame()

for(cont in c('africa','antarctica','asia','europe','northamerica','oceania','southamerica')){
  cont_dir <- paste0('output/ramsar_metric_results_100m/', cont)
  data_dir <- paste0('inst/extdata/ramsar_site_data_100m_', cont)
  wkt_dir <- paste0('inst/extdata/ramsar_sites_wkt')
  
  if(!dir.exists(cont_dir)) next
  
  files <- list.files(cont_dir, pattern='total_occ_map.RData$', recursive=TRUE)
  
  for(f in files){
    obj <- tryCatch({
      readRDS(paste0(cont_dir, '/', f))
    }, error = function(e) NULL)
    
    if(is.null(obj) || is.null(obj$data)) next
    
    # Extract site_id and country from filename
    site_id <- sub('_total_occ_map.RData$', '', basename(f))
    country <- dirname(f)
    country <- basename(country)
    
    # Load the CSV to filter by year
    csv_path <- file.path(data_dir, country, paste0(site_id, '_data.csv'))
    
    if(!file.exists(csv_path)) next
    
    tryCatch({
      dat <- read.csv(csv_path)
      
      # Filter to 1990 onwards
      dat_1990 <- dat[dat$year >= year_threshold, ]
      occ_1990 <- sum(dat_1990$occurrences, na.rm = TRUE)
      
      # Get total for comparison
      occ_all <- sum(dat$occurrences, na.rm = TRUE)
      
    }, error = function(e) {
      occ_1990 <<- NA
      occ_all <<- NA
    })
    
    if(is.na(occ_1990) || occ_1990 == 0) next
    
    # Get area from WKT file
    wkt_file <- file.path(wkt_dir, country, paste0(site_id, '.wkt'))
    
    actual_area <- NA
    
    if(file.exists(wkt_file)){
      tryCatch({
        wkt <- paste(readLines(wkt_file), collapse='')
        poly <- st_as_sfc(wkt)
        st_crs(poly) <- 4326
        actual_area <- as.numeric(st_area(poly)) / 1e6
      }, error = function(e) NA)
    }
    
    # Fall back to bbox
    bbox_area <- NA
    if(is.na(actual_area)){
      tryCatch({
        bbox_area <- as.numeric(st_area(obj$original_bbox)) / 1e6
      }, error = function(e) NA)
    }
    
    area <- if(!is.na(actual_area)) actual_area else bbox_area
    
    if(is.na(area) || area <= 0) next
    
    # Calculate density
    density_1990 <- occ_1990 / area
    density_all <- occ_all / area
    
    results <- rbind(results, data.frame(
      continent=cont,
      site_id=site_id,
      country=country,
      density_1990=density_1990,
      density_all=density_all,
      passes_1990=as.numeric(density_1990 >= 0.25),
      passes_all=as.numeric(density_all >= 0.25),
      area_km2=area,
      occ_1990=occ_1990,
      occ_all=occ_all,
      area_source=ifelse(!is.na(actual_area), "WKT", "bbox"),
      stringsAsFactors=FALSE
    ))
  }
  
  cat('Processed', cont, ':', nrow(results), 'sites total\n')
}

# Save full results
write.csv(results, 'output/data_sufficiency/density_analysis_1990.csv', row.names=FALSE)

total <- nrow(results)
passing_1990 <- sum(results$passes_1990 == 1, na.rm=TRUE)
passing_all <- sum(results$passes_all == 1, na.rm=TRUE)

cat('\n===== DATA SUFFICIENCY ANALYSIS (>=0.25/km2) =====\n')
cat('Total sites:', total, '\n')
cat('Sites passing (1990+):', passing_1990, '(', round(passing_1990/total*100,1), '%)\n')
cat('Sites passing (all time):', passing_all, '(', round(passing_all/total*100,1), '%)\n\n')

cat('By continent (1990+):\n')
for (cont in c('africa','antarctica','asia','europe','northamerica','oceania','southamerica')) {
  sub <- results[results$continent==cont,]
  n <- nrow(sub)
  p <- sum(sub$passes_1990 == 1, na.rm=TRUE)
  if(n>0) cat(sprintf('  %s: %d sites, %d passing (%.1f%%)\n', cont, n, p, p/n*100))
}

# Example
cat('\nExample - Algeria 1052:\n')
print(subset(results, grepl('1052', site_id) & continent == 'africa'))

# Save passing (1990+)
passing <- results[results$passes_1990 == 1,]
write.csv(passing, 'output/data_sufficiency/sites_for_b3gbi_1990.csv', row.names=FALSE)

cat('\nSaved to output/data_sufficiency/\n')
