#!/usr/bin/env Rscript
# =============================================================================
# run_temporal_only.R
#
# Drives worker_temporal_only.R over every site cube.
#
# NON-DESTRUCTIVE: writes only output/data_sufficiency/subgroup_analysis/
# temporal_only_results.csv and uses temp_workers_temporal/ as scratch.
# Reads nothing that it writes. Resumable.
#
# Site discovery and de-duplication are copied verbatim from
# run_subgroup_data_sufficiency_parallel.R so the same set of cubes is used.
# (Checked: no site_id occurs under more than one continent/country, so the
# group_by(site_id) de-duplication cannot cross-match.)
# =============================================================================
library(dplyr)
library(parallel)

CONTINENTS <- c("africa", "antarctica", "asia", "europe", "northamerica", "oceania", "southamerica")
OUT_NAME <- "temporal_only_results.csv"

run_master <- function() {
  out_dir <- "output/data_sufficiency/subgroup_analysis"
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_file <- file.path(out_dir, OUT_NAME)

  if (grepl("raw_tests_results\\.csv$|all_tests_results", out_file)) {
    stop("Refusing to write to an original results file.")
  }

  processed_sites <- c()
  if (file.exists(out_file)) {
    prev <- tryCatch(read.csv(out_file, stringsAsFactors = FALSE), error = function(e) data.frame())
    if (nrow(prev) > 0 && "site_id" %in% names(prev)) processed_sites <- unique(prev$site_id)
  } else {
    cat("continent,country,site_id,subset,temporal_distance_cumsum,temporal_distance_annual,total_occurrences\n",
        file = out_file)
  }
  cat(sprintf("Already processed: %d sites.\n", length(processed_sites)))

  all_files_df <- data.frame(cont = character(), country = character(), site_id = character(),
                             file_path = character(), size = numeric(), stringsAsFactors = FALSE)

  for (cont in CONTINENTS) {
    data_dir <- file.path("inst/extdata", paste0("ramsar_site_data_100m_", cont))
    if (!dir.exists(data_dir)) next
    countries <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
    for (country in countries) {
      cube_dir <- file.path(data_dir, country)
      for (cf in list.files(cube_dir, pattern = "\\.csv$", full.names = TRUE)) {
        fn <- basename(cf)
        m <- regexpr("^site_[0-9]+", fn)
        site_id <- if (m != -1) regmatches(fn, m) else sub("(_[a-zA-Z0-9_]+)?_data\\.csv$", "", fn)
        all_files_df <- rbind(all_files_df, data.frame(
          cont = cont, country = country, site_id = site_id,
          file_path = cf, size = file.info(cf)$size, stringsAsFactors = FALSE))
      }
    }
  }

  deduped <- all_files_df %>%
    group_by(site_id) %>% filter(size == max(size)) %>% filter(row_number() == 1) %>% ungroup()

  site_list <- list()
  for (i in seq_len(nrow(deduped))) {
    row <- deduped[i, ]
    if (!(row$site_id %in% processed_sites)) {
      site_list[[length(site_list) + 1]] <- list(cont = row$cont, country = row$country, site_id = row$site_id)
    }
  }

  total_sites <- length(site_list)
  cat(sprintf("Sites remaining: %d\n", total_sites))
  if (total_sites == 0) { cat("Nothing to do.\n"); return() }

  num_cores <- 8
  temp_dir <- file.path(out_dir, "temp_workers_temporal")
  if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
  dir.create(temp_dir, recursive = TRUE)

  cl <- makePSOCKcluster(num_cores)
  clusterExport(cl, c("temp_dir", "site_list"), envir = environment())

  results <- parLapply(cl, seq_along(site_list), function(i) {
    site <- site_list[[i]]
    worker_out <- file.path(temp_dir, paste0("worker_", i, ".csv"))
    res <- system2("Rscript", args = c("scripts/worker_temporal_only.R",
                                       shQuote(site$cont), shQuote(site$country),
                                       shQuote(site$site_id), shQuote(worker_out)),
                   stdout = FALSE, stderr = FALSE)
    return(res == 0)
  })
  stopCluster(cl)

  cat(sprintf("Finished. Failures: %d\n", sum(unlist(results) == FALSE)))

  merged <- 0
  for (tf in list.files(temp_dir, pattern = "worker_.*\\.csv$", full.names = TRUE)) {
    lines <- readLines(tf, warn = FALSE)
    if (length(lines) > 0) { cat(lines, sep = "\n", file = out_file, append = TRUE); merged <- merged + 1 }
  }
  unlink(temp_dir, recursive = TRUE)

  cat(sprintf("Merged %d worker files into %s\n", merged, OUT_NAME))
  cat("Next: Rscript scripts/merge_and_regate_v2.R\n")
}

run_master()
