# ============================================================================
# SCRIPT: Robust Automated Continental GBIF Data Cube Download
#
# !!! CRITICAL: MASSIVE DATA VOLUME WARNING !!!
# Continental GBIF downloads (e.g., ASIA, EUROPE) can exceed 100GB.
#
# REQUIREMENTS:
# 1. Disk Space: At least 200GB free (for compressed download + extraction).
# 2. Internet: High-speed, stable connection (resuming is supported).
# 3. Memory: This script is low-memory, but 'scripts/split_data_cubes.R' 
#    requires significant RAM or small chunk sizes.
# ============================================================================

# --- 0. Setup ---
lapply(list.files("R", full.names = TRUE), source)
library(rgbif)
library(dplyr)

# --- Configuration ---
target_continent <- "ASIA"  # 'AFRICA', 'ANTARCTICA', 'ASIA', 'EUROPE', etc.
mgrs_resolution <- 100      # 100m resolution
output_dir <- file.path("inst/extdata/continental_gbif_data", tolower(target_continent))
key_file <- file.path(output_dir, ".download_key")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# --- Pre-flight Check ---
check_gbif_credentials()

# --- 1. Recovery Check ---
# Check if we have an existing job running
existing_key <- load_download_key(key_file)
download_key <- NULL

if (!is.null(existing_key)) {
  message("Found existing download key: ", existing_key)
  meta <- occ_download_meta(existing_key)
  
  if (meta$status %in% c("RUNNING", "PREPARING", "SUCCEEDED")) {
    user_resume <- readline(prompt = paste0("Job '", existing_key, "' is ", meta$status, ". Resume this job? (yes/no): "))
    if (tolower(user_resume) == "yes") {
      download_key <- existing_key
    }
  }
}

# --- 2. Request New Job (if not resuming) ---
if (is.null(download_key)) {
  sql_query <- sprintf(
    "\n    SELECT\n      GBIF_MGRSCode(%d, decimalLatitude, decimalLongitude, 0.) AS mgrscellcode,\n      year,\n      speciesKey,\n      datasetKey,\n      COUNT(*) AS occurrences\n    FROM\n      occurrence\n    WHERE\n      continent = '%s'\n      AND hasCoordinate = TRUE\n      AND speciesKey IS NOT NULL\n      AND occurrenceStatus = 'PRESENT'\n    GROUP BY\n      mgrscellcode,\n      year,\n      speciesKey,\n      datasetKey\n  ", mgrs_resolution, target_continent)

  message("\n--- NEW CONTINENTAL REQUEST ---")
  message("Continent: ", target_continent, " | Resolution: ", mgrs_resolution, "m")
  message("EXPECTED SIZE: >100GB")
  
  confirm <- readline(prompt = "Submit this massive SQL job to GBIF? (yes/no): ")
  if (tolower(confirm) != "yes") stop("Cancelled by user.")
  
  message("Submitting SQL request...")
  res <- occ_download_sql(q = sql_query)
  download_key <- res[1]
  save_download_key(download_key, key_file)
}

# --- 3. Wait for GBIF Processing ---
message("\nMonitoring job status (Link: https://www.gbif.org/occurrence/download/", download_key, ")")
status <- occ_download_wait(download_key)

if (status$status != "SUCCEEDED") {
  stop("GBIF job failed or was cancelled. Status: ", status$status)
}

# --- 4. Resumable File Retrieval ---
# We use the direct download URL + system curl for resume support
download_url <- paste0("https://api.gbif.org/v1/occurrence/download/request/", download_key, ".zip")
dest_zip <- file.path(output_dir, paste0(tolower(target_continent), "_cube.zip"))

message("\n--- STARTING FILE RETRIEVAL ---")
message("If the download cuts out, simply restart this script to resume.")

# Robust download via curl
download_robust(download_url, dest_zip)

# --- 5. Extraction ---
# Note: unzip() in R may fail for files >4GB on some systems. 
# If it fails, use a system tool like 7-Zip or tar.
message("\n--- EXTRACTION ---")
message("Extracting massive archive... this will take time.")
tryCatch({
  utils::unzip(dest_zip, exdir = output_dir)
  message("Extraction complete. You can now delete the .zip to save space.")
}, error = function(e) {
  message("R's unzip failed (common for >4GB files). Please manually extract: ", dest_zip)
  message("Error: ", e$message)
})

message("\nWorkflow Finished for: ", target_continent)
