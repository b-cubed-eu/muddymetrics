#!/usr/bin/env Rscript

results <- data.frame()

for(cont in c('africa','antarctica','asia','europe','northamerica','oceania','southamerica')){
  cont_dir <- paste0('output/ramsar_metric_results_100m/', cont)
  if(!dir.exists(cont_dir)) next
  
  # Find all overall_density.RData files recursively
  files <- list.files(cont_dir, pattern='overall_density.RData$', recursive=TRUE)
  
  for(f in files){
    # Read with readRDS
    obj <- tryCatch({
      readRDS(paste0(cont_dir, '/', f))
    }, error = function(e) NULL)
    
    if(is.null(obj) || length(obj) < 1 || is.null(obj[[1]])) next
    
    # Extract density value from units object
    density_val <- tryCatch(as.numeric(obj[[1]]), error=function(e) NA)
    if(is.na(density_val)) next
    
    # Extract site_id from filename
    site_id <- sub('_overall_density.RData$', '', basename(f))
    
    # Extract country from path
    country <- dirname(f)
    country <- basename(country)
    
    results <- rbind(results, data.frame(
      continent=cont,
      site_id=site_id,
      country=country,
      density=density_val,
      passes=as.numeric(density_val >= 0.25),
      stringsAsFactors=FALSE
    ))
  }
  
  cat('Processed', cont, ':', nrow(results), 'sites total\n')
}

# Save
write.csv(results, 'output/data_sufficiency/density_analysis.csv', row.names=FALSE)

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
print(subset(results, grepl('1052', site_id)))

# Save passing
passing <- results[results$passes == 1,]
write.csv(passing, 'output/data_sufficiency/sites_for_b3gbi.csv', row.names=FALSE)

cat('\nSaved to output/data_sufficiency/\n')
