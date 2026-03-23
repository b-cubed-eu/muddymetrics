#!/usr/bin/env Rscript
library(sf)
library(units)

results <- data.frame()

for(cont in c('africa','antarctica','asia','europe','northamerica','oceania','southamerica')){
  cont_dir <- paste0('output/ramsar_metric_results_100m/', cont)
  if(!dir.exists(cont_dir)) next
  
  files <- list.files(cont_dir, pattern='total_occ_map.RData$', recursive=TRUE)
  
  for(f in files){
    obj <- tryCatch({
      readRDS(paste0(cont_dir, '/', f))
    }, error = function(e) NULL)
    
    if(is.null(obj) || is.null(obj$original_bbox) || is.null(obj$data)) next
    
    # Calculate area
    area <- tryCatch({
      sf::st_area(obj$original_bbox) |> set_units("km^2") |> as.numeric()
    }, error = function(e) NA)
    
    if(is.na(area) || area <= 0) next
    
    # Calculate occurrence sum
    occ <- sum(obj$data$diversity_val, na.rm = TRUE)
    
    # Calculate density
    density_val <- occ / area
    
    # Extract site_id and country from filename
    site_id <- sub('_total_occ_map.RData$', '', basename(f))
    country <- dirname(f)
    country <- basename(country)
    
    results <- rbind(results, data.frame(
      continent=cont,
      site_id=site_id,
      country=country,
      density=density_val,
      passes=as.numeric(density_val >= 0.25),
      area_km2=area,
      occ_count=occ,
      stringsAsFactors=FALSE
    ))
  }
  
  cat('Processed', cont, ':', nrow(results), 'sites total\n')
}

# Save full results
write.csv(results, 'output/data_sufficiency/density_analysis_corrected.csv', row.names=FALSE)

total <- nrow(results)
passing_count <- sum(results$passes == 1, na.rm=TRUE)

cat('\n===== DATA SUFFICIENCY ANALYSIS =====\n')
cat('Total sites:', total, '\n')
cat('Sites passing (>= 0.25/km2):', passing_count, '\n')
cat('Pass rate:', round(passing_count/total*100,1), '%\n\n')

cat('By continent:\n')
for (cont in c('africa','antarctica','asia','europe','northamerica','oceania','southamerica')) {
  sub <- results[results$continent==cont,]
  n <- nrow(sub)
  p <- sum(sub$passes == 1, na.rm=TRUE)
  if(n>0) cat(sprintf('  %s: %d sites, %d passing (%.1f%%)\n', cont, n, p, p/n*100))
}

# Example
cat('\nExample - Algeria 1052:\n')
print(subset(results, grepl('1052', site_id) & continent == 'africa'))

# Save passing
passing <- results[results$passes == 1,]
write.csv(passing, 'output/data_sufficiency/sites_for_b3gbi_corrected.csv', row.names=FALSE)

cat('\nSaved to output/data_sufficiency/\n')
