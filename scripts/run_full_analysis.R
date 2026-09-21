#!/usr/bin/env Rscript

library(callr)
library(dplyr)

out_dir <- "output/full_analysis"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

in_dir <- "output/data_sufficiency/subgroup_analysis"

file_passing <- file.path(in_dir, "sites_passing_corrected.csv")
job_queue <- list()

if (file.exists(file_passing)) {
  df <- read.csv(file_passing, stringsAsFactors = FALSE)
  if (nrow(df) > 0) {
    for (i in 1:nrow(df)) {
      job_queue[[i]] <- list(
        cont = df$continent[i],
        country = df$country[i],
        site_id = df$site_id[i],
        subset = df$subset[i],
        run_ts = TRUE,
        run_map = FALSE
      )
    }
  }
}

total_jobs <- length(job_queue)
cat(sprintf("Found %d site-subgroup combinations to process for full analysis.\n", total_jobs))
flush.console()

# Clear the progress log at start
writeLines(character(0), "full_analysis_progress.log")


library(parallel)

run_job <- function(job) {
  # Log progress to a shared file
  cat(sprintf("[%s] Starting: %s -> %s (%s)\n", format(Sys.time(), "%H:%M:%S"), job$country, job$site_id, job$subset),
      file = "full_analysis_progress.log", append = TRUE)
  
  out_dir <- "output/full_analysis"
  site_dir <- file.path(out_dir, job$site_id)
  if (!dir.exists(site_dir)) dir.create(site_dir, recursive = TRUE)
  log_file <- file.path(site_dir, paste0("worker_", job$subset, ".log"))
  
  res <- tryCatch({
    callr::rscript("scripts/worker_full_analysis.R", 
                   cmdargs = c(job$cont, job$country, job$site_id, job$subset, 
                               as.character(job$run_ts), as.character(job$run_map)),
                   timeout = 1800,
                   stdout = log_file,
                   stderr = log_file,
                   show = FALSE)
    0  # success
  }, callr_timeout_error = function(e) {
    99 # timeout
  }, error = function(e) {
    1  # error
  })
  
  cat(sprintf("[%s] Finished (code %d): %s -> %s (%s)\n", format(Sys.time(), "%H:%M:%S"), res, job$country, job$site_id, job$subset),
      file = "full_analysis_progress.log", append = TRUE)
  return(res)
}

num_cores <- min(detectCores() - 1, 2)
cat(sprintf("Running in parallel using %d cores. Monitoring logs in 'full_analysis_progress.log'...\n", num_cores))
flush.console()

cl <- makePSOCKcluster(num_cores)
clusterExport(cl, c("run_job"))
clusterEvalQ(cl, { library(callr) })

results_list <- parLapply(cl, job_queue, run_job)
stopCluster(cl)

results <- unlist(results_list)
failures <- sum(results == 1)
timeouts <- sum(results == 99)

cat(sprintf("Full analysis completed. Processed %d jobs. Crashes: %d, Timeouts: %d.\n", total_jobs, failures, timeouts))
