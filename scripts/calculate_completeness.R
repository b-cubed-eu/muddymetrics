#!/usr/bin/env Rscript
library(b3gbi)
library(readr)

passing_sites <- read.csv('output/data_sufficiency/standard_0.25/sites_passing.csv')

cat("Total passing sites:", nrow(passing_sites), "\n\n")

results <- data.frame()

for(cont in unique(passing_sites$continent)){
  cont_sites <- passing_sites[passing_sites$continent == cont,]
  data_dir <- paste0("inst/extdata/ramsar_site_data_100m_", cont)
  
  if(!dir.exists(data_dir)){
    cat("Skipping", cont, "\n")
    next
  }
  
  cat("Processing continent:", cont, "-", nrow(cont_sites), "sites\n")
  
  for(idx in 1:nrow(cont_sites)){
    site_id <- cont_sites$site_id[idx]
    country <- cont_sites$country[idx]
    
    cubepath <- file.path(data_dir, country, paste0(site_id, "_data.csv"))
    
    if(!file.exists(cubepath)){
      next
    }
    
    tryCatch({
      cube <- b3gbi::process_cube(cubepath)
      result <- b3gbi::completeness_ts(cube, ci_type='none')
      
      latest <- tail(result$data$diversity_val, 1)
      mean_val <- mean(result$data$diversity_val, na.rm=TRUE)
      
      results <<- rbind(results, data.frame(
        continent = cont,
        site_id = site_id,
        country = country,
        completeness_latest = latest,
        completeness_mean = mean_val,
        passes = as.numeric(!is.na(latest) && latest >= 0.70),
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e){
      # Skip on error
    })
    
    if(idx %% 50 == 0) cat("  Processed", idx, "of", nrow(cont_sites), "\n")
  }
}

# Save
write.csv(results, 'output/data_sufficiency/chao2_sac/completeness_analysis.csv', row.names=FALSE)

# Summary
cat('\n===== COMPLETENESS ANALYSIS =====\n')
cat('Sites with results:', nrow(results), '\n')
cat('Passing (>=70%):', sum(results$passes == 1), '\n\n')

cat('By continent:\n')
for(cont in c('africa','asia','europe','northamerica','oceania','southamerica')){
  sub <- results[results$continent == cont,]
  n <- nrow(sub)
  p <- sum(sub$passes == 1, na.rm=TRUE)
  if(n > 0) cat(sprintf('  %s: %d sites, %d passing (%.1f%%)\n', cont, n, p, p/n*100))
}

cat('\nCompleteness distribution:\n')
cat('  Mean:', round(mean(results$completeness_latest, na.rm=TRUE), 3), '\n')
cat('  Median:', round(median(results$completeness_latest, na.rm=TRUE), 3), '\n')
"