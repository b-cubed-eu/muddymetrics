#!/usr/bin/env Rscript
# Calculate completeness using b3gbi::completeness_ts() - parallel version

library(b3gbi)
library(parallel)
library(foreach)
library(doParallel)

# Read passing sites
passing_sites <- read.csv('output/data_sufficiency/standard_0.25/sites_passing.csv')

cat("Total passing sites:", nrow(passing_sites), "\n")

# Setup parallel
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
registerDoParallel(cl)

threshold_completeness <- 0.70

# Create work items
work_items <- data.frame()

for(cont in unique(passing_sites$continent)){
  cont_sites <- passing_sites[passing_sites$continent == cont,]
  data_dir <- paste0("inst/extdata/ramsar_site_data_100m_", cont)
  wkt_dir <- "inst/extdata/ramsar_sites_wkt"
  
  if(!dir.exists(data_dir)) next
  
  for(idx in 1:nrow(cont_sites)){
    work_items <- rbind(work_items, data.frame(
      continent = cont,
      site_id = cont_sites$site_id[idx],
      country = cont_sites$country[idx],
      data_dir = data_dir,
      wkt_dir = wkt_dir
    ))
  }
}

cat("Total work items:", nrow(work_items), "\n")

# Process in parallel
results <- foreach(i = 1:nrow(work_items), .combine = rbind, .packages = c("b3gbi")) %dopar% {
  item <- work_items[i,]
  
  cubepath <- file.path(item$data_dir, item$country, paste0(item$site_id, "_data.csv"))
  shapefilepath <- file.path(item$wkt_dir, item$country, paste0(item$site_id, ".wkt"))
  
  if(!file.exists(cubepath) || !file.exists(shapefilepath)){
    return(NULL)
  }
  
  tryCatch({
    cube <- b3gbi::process_cube(cubepath, separator = ",")
    
    comp_result <- b3gbi::compute_indicator_workflow(
      data = cube,
      type = "completeness",
      dim_type = "ts",
      shapefile_path = shapefilepath,
      ne_scale = "large",
      region = item$continent,
      include_water = TRUE,
      shapefile_crs = 4326,
      ci_type = "none"
    )
    
    comp_values <- comp_result$data$diversity_val
    latest <- tail(comp_values[!is.na(comp_values)], 1)
    if(length(latest) == 0) latest <- NA
    
    data.frame(
      continent = item$continent,
      site_id = item$site_id,
      country = item$country,
      completeness_latest = latest,
      passes = as.numeric(!is.na(latest) && latest >= 0.70),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

stopCluster(cl)

# Save
write.csv(results, 'output/data_sufficiency/chao2_sac/completeness_analysis.csv', row.names = FALSE)

# Summary
cat('\n===== COMPLETENESS ANALYSIS =====\n')
cat('Sites with results:', nrow(results), '\n')
cat('Passing (>=70%):', sum(results$passes == 1), '\n')

cat('\nBy continent:\n')
for(cont in c('africa','asia','europe','northamerica','oceania','southamerica')){
  sub <- results[results$continent == cont,]
  n <- nrow(sub)
  p <- sum(sub$passes == 1, na.rm = TRUE)
  if(n > 0) cat(sprintf('  %s: %d sites, %d passing (%.1f%%)\n', cont, n, p, p/n*100))
}
