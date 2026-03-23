#!/usr/bin/env Rscript
# Calculate Chao2 and SAC slope for sites passing data sufficiency
# SAC slope based on temporal accumulation (years)

library(readr)

calculate_chao2 <- function(cube) {
  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }
  
  species_col <- cube$specieskey
  valid_species <- species_col[!is.na(species_col)]
  
  if (length(valid_species) == 0) {
    return(list(
      chao2 = NA_real_,
      observed = 0,
      completeness = NA_real_,
      f1 = 0,
      f2 = 0
    ))
  }
  
  sample_counts <- table(valid_species)
  
  f1 <- sum(sample_counts == 1)
  f2 <- sum(sample_counts == 2)
  
  S_obs <- length(unique(valid_species))
  
  if (f2 > 0) {
    chao2_est <- S_obs + (f1^2) / (2 * f2)
  } else if (f1 > 0) {
    chao2_est <- S_obs + (f1 * (f1 - 1)) / 2
  } else {
    chao2_est <- S_obs
  }
  
  completeness <- if (chao2_est > 0) S_obs / chao2_est else NA_real_
  
  return(list(
    chao2 = chao2_est,
    observed = S_obs,
    completeness = completeness,
    f1 = f1,
    f2 = f2
  ))
}

calculate_sac_slope <- function(cube, n_iterations = 50) {
  if (is.character(cube) && file.exists(cube)) {
    cube <- readr::read_csv(cube, show_col_types = FALSE)
  }
  
  # Get unique years
  years <- sort(unique(cube$year))
  
  if (length(years) < 3) {
    return(list(
      slope = NA_real_,
      intercept = NA_real_,
      r_squared = NA_real_
    ))
  }
  
  tryCatch({
    # Calculate cumulative species over time
    species_accum <- numeric()
    year_accum <- numeric()
    
    for (i in 1:length(years)) {
      sub_cube <- cube[cube$year <= years[i], ]
      n_species <- length(unique(sub_cube$specieskey[!is.na(sub_cube$specieskey) & sub_cube$specieskey > 0]))
      species_accum <- c(species_accum, n_species)
      year_accum <- c(year_accum, years[i])
    }
    
    # Normalize years to 0-1 scale for slope interpretation
    year_norm <- (year_accum - min(year_accum)) / (max(year_accum) - min(year_accum) + 1)
    
    if (length(species_accum) < 3) {
      return(list(slope = NA_real_, intercept = NA_real_, r_squared = NA_real_))
    }
    
    lm_result <- lm(species_accum ~ year_norm)
    r_squared <- summary(lm_result)$r.squared
    
    # The slope represents species accumulation rate per unit time
    # Lower slope = species accumulation slowing = more complete sampling
    slope <- coef(lm_result)[2]
    
    return(list(
      slope = slope,
      intercept = coef(lm_result)[1],
      r_squared = r_squared
    ))
  }, error = function(e) {
    return(list(slope = NA_real_, intercept = NA_real_, r_squared = NA_real_))
  })
}

# Read passing sites (from class-based density)
passing_sites <- read.csv('output/data_sufficiency/sites_for_b3gbi_by_class_1990.csv')

cat("Total passing sites:", nrow(passing_sites), "\n")

# Define thresholds from Troia 2016
threshold_chao2 <- 0.70
threshold_slope <- 0.10

results <- data.frame()

for(cont in unique(passing_sites$continent)){
  cont_sites <- passing_sites[passing_sites$continent == cont,]
  data_dir <- paste0("inst/extdata/ramsar_site_data_100m_", cont)
  
  if(!dir.exists(data_dir)){
    cat("Skipping", cont, "- no data dir\n")
    next
  }
  
  cat("Processing continent:", cont, "\n")
  
  for(idx in 1:nrow(cont_sites)){
    site_row <- cont_sites[idx,]
    site_id <- site_row$site_id
    country <- site_row$country
    
    cubepath <- file.path(data_dir, country, paste0(site_id, "_data.csv"))
    
    if(!file.exists(cubepath)){
      cat("  Missing cube:", site_id, "\n")
      next
    }
    
    tryCatch({
      # Filter to 1990+
      dat <- read.csv(cubepath)
      dat_1990 <- dat[dat$year >= 1990, ]
      
      # Skip if too few records
      if(nrow(dat_1990) < 10) next
      
      # Write temp filtered cube
      temp_cube <- tempfile(fileext = ".csv")
      write.csv(dat_1990, temp_cube, row.names = FALSE)
      
      # Calculate chao2
      chao2_result <- calculate_chao2(temp_cube)
      
      # Calculate SAC slope
      sac_result <- calculate_sac_slope(temp_cube)
      
      # Check thresholds
      passes_chao2 <- !is.na(chao2_result$completeness) && chao2_result$completeness >= threshold_chao2
      passes_slope <- !is.na(sac_result$slope) && sac_result$slope <= threshold_slope
      
      results <- rbind(results, data.frame(
        continent = cont,
        site_id = site_id,
        country = country,
        chao2 = chao2_result$chao2,
        observed_species = chao2_result$observed,
        completeness_chao2 = chao2_result$completeness,
        passes_chao2 = as.numeric(passes_chao2),
        sac_slope = sac_result$slope,
        sac_r2 = sac_result$r_squared,
        passes_slope = as.numeric(passes_slope),
        passes_both = as.numeric(passes_chao2 & passes_slope),
        stringsAsFactors = FALSE
      ))
      
      unlink(temp_cube)
      
    }, error = function(e){
      cat("  Error:", site_id, "-", conditionMessage(e), "\n")
    })
    
    if(idx %% 100 == 0) cat("  Processed", idx, "of", nrow(cont_sites), "\n")
  }
}

# Save results
write.csv(results, 'output/data_sufficiency/chao2_sac_analysis.csv', row.names = FALSE)

# Summary
cat('\n===== CHAO2 & SAC SLOPE ANALYSIS =====\n')
cat('Total sites:', nrow(results), '\n')
cat('Passing Chao2 (>=70%):', sum(results$passes_chao2 == 1), '\n')
cat('Passing SAC slope (<=0.10):', sum(results$passes_slope == 1), '\n')
cat('Passing BOTH:', sum(results$passes_both == 1), '\n')

cat('\nSAC slope distribution:\n')
cat('  Mean:', mean(results$sac_slope, na.rm=TRUE), '\n')
cat('  Median:', median(results$sac_slope, na.rm=TRUE), '\n')
cat('  Min:', min(results$sac_slope, na.rm=TRUE), '\n')
cat('  Max:', max(results$sac_slope, na.rm=TRUE), '\n')

cat('\nBy continent:\n')
for(cont in c('africa','antarctica','asia','europe','northamerica','oceania','southamerica')){
  sub <- results[results$continent == cont,]
  n <- nrow(sub)
  p_chao <- sum(sub$passes_chao2 == 1, na.rm = TRUE)
  p_slope <- sum(sub$passes_slope == 1, na.rm = TRUE)
  p_both <- sum(sub$passes_both == 1, na.rm = TRUE)
  if(n > 0) cat(sprintf('  %s: %d sites | Chao2: %d | Slope: %d | Both: %d\n', 
                        cont, n, p_chao, p_slope, p_both))
}

# Example
cat('\nExample - some Europe sites:\n')
print(head(subset(results, continent == 'europe'), 10))

# Save sites passing both thresholds
passing_both <- results[results$passes_both == 1, ]
write.csv(passing_both, 'output/data_sufficiency/sites_chao2_sac_pass.csv', row.names = FALSE)

cat('\nSaved to output/data_sufficiency/\n')
