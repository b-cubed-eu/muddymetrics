#!/usr/bin/env Rscript
library(dplyr)
library(parallel)

CONTINENTS <- c("africa", "antarctica", "asia", "europe", "northamerica", "oceania", "southamerica")

run_master_parallel <- function() {
  out_dir <- "output/data_sufficiency/subgroup_analysis"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  out_file <- file.path(out_dir, "raw_tests_results.csv")
  
  # Check already processed sites in out_file
  processed_sites <- c()
  if (file.exists(out_file)) {
    raw_res <- tryCatch(read.csv(out_file, stringsAsFactors = FALSE), error = function(e) data.frame())
    if (nrow(raw_res) > 0 && "site_id" %in% names(raw_res)) {
      processed_sites <- unique(raw_res$site_id)
    }
  } else {
    cat("continent,country,site_id,subset,density,chao2,sac_slope,temporal_distance,spatial_distance\n", file = out_file)
  }
  
  cat(sprintf("Found %d already processed sites in raw results.\n", length(processed_sites)))
  
  # Gather all sites from directories and deduplicate based on file size
  all_files_df <- data.frame(
    cont = character(),
    country = character(),
    site_id = character(),
    file_path = character(),
    size = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cont in CONTINENTS) {
    data_dir <- file.path("inst/extdata", paste0("ramsar_site_data_100m_", cont))
    if (!dir.exists(data_dir)) next
    
    countries <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
    for (country in countries) {
      cube_dir <- file.path(data_dir, country)
      cube_files <- list.files(cube_dir, pattern = "\\.csv$", full.names = TRUE)
      
      for (cf in cube_files) {
        fn <- basename(cf)
        site_id_match <- regexpr("^site_[0-9]+", fn)
        if (site_id_match != -1) {
            site_id <- regmatches(fn, site_id_match)
        } else {
            site_id <- sub("(_[a-zA-Z0-9_]+)?_data\\.csv$", "", fn)
        }
        
        sz <- file.info(cf)$size
        all_files_df <- rbind(all_files_df, data.frame(
          cont = cont,
          country = country,
          site_id = site_id,
          file_path = cf,
          size = sz,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  # Keep only the row with the largest file size for each unique site_id
  deduped_files <- all_files_df %>%
    group_by(site_id) %>%
    filter(size == max(size)) %>%
    filter(row_number() == 1) %>%
    ungroup()
  
  site_list <- list()
  for (i in seq_len(nrow(deduped_files))) {
    row <- deduped_files[i, ]
    if (!(row$site_id %in% processed_sites)) {
      site_list[[length(site_list) + 1]] <- list(
        cont = row$cont,
        country = row$country,
        site_id = row$site_id
      )
    }
  }
  
  total_sites <- length(site_list)
  cat(sprintf("Found %d sites remaining to process.\n", total_sites))
  
  if (total_sites == 0) {
    cat("All sites are already processed!\n")
    return()
  }
  
  # Setup parallel cluster
  num_cores <- 8
  cat(sprintf("Starting parallel execution using %d cores...\n", num_cores))
  
  temp_dir <- file.path(out_dir, "temp_workers")
  if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
  dir.create(temp_dir, recursive = TRUE)
  
  cl <- makePSOCKcluster(num_cores)
  clusterExport(cl, c("temp_dir", "site_list"), envir = environment())
  
  # Run workers
  results <- parLapply(cl, seq_along(site_list), function(i) {
    site <- site_list[[i]]
    worker_out <- file.path(temp_dir, paste0("worker_", i, ".csv"))
    
    # We pass worker_out as the destination file for this single site
    res <- system2("Rscript", args = c("scripts/worker_subgroup.R", 
                                       shQuote(site$cont), 
                                       shQuote(site$country), 
                                       shQuote(site$site_id), 
                                       shQuote(worker_out)), 
                   stdout = FALSE, stderr = FALSE)
    return(res == 0)
  })
  
  stopCluster(cl)
  
  failures <- sum(unlist(results) == FALSE)
  cat(sprintf("Finished parallel execution. Failures: %d\n", failures))
  
  # Merge individual site results
  cat("Merging new results into raw_tests_results.csv...\n")
  temp_files <- list.files(temp_dir, pattern = "worker_.*\\.csv$", full.names = TRUE)
  
  merged_count <- 0
  for (tf in temp_files) {
    lines <- readLines(tf, warn = FALSE)
    if (length(lines) > 0) {
      cat(lines, sep = "\n", file = out_file, append = TRUE)
      merged_count <- merged_count + 1
    }
  }
  
  # Clean up temp files
  unlink(temp_dir, recursive = TRUE)
  
  cat(sprintf("Successfully merged %d files. Running compute_thresholds.R to update passing lists...\n", merged_count))
  system2("Rscript", args = c("scripts/compute_thresholds.R"))
  cat("Master parallel run completed successfully.\n")
}

run_master_parallel()
