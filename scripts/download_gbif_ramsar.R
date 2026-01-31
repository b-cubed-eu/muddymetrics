# --- Setup ---
lapply(list.files("R", full.names = TRUE), source)
library(rgbif)
library(sf)
library(dplyr)
library(utils)

# --- Pre-flight Check ---
check_gbif_credentials()

# --- Configuration ---
# Base directory where your country-specific Ramsar site folders (containing WKT files) are.
base_ramsar_dir <- "inst/extdata/ramsar_sites_wkt"

# Directory to save the downloaded GBIF data
output_gbif_base_dir <- "inst/extdata/gbif_ramsar" 
if (!dir.exists(output_gbif_base_dir)) {
  dir.create(output_gbif_base_dir, recursive = TRUE)
}

# --- GBIF API parameters ---
gbif_base_filters <- list(
  hasGeospatialIssue = FALSE,
  hasCoordinate = TRUE,
  basisOfRecord = c("OBSERVATION", "HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "LITERATURE", "LIVING_SPECIMEN", "MATERIAL_SAMPLE", "MACHINE_OBSERVATION")
)

# --- Get all country subdirectories ---
country_dirs <- list.dirs(base_ramsar_dir, recursive = FALSE, full.names = TRUE)

if (length(country_dirs) == 0) {
  stop("No country subdirectories found in: ", base_ramsar_dir)
}

message(paste("Found", length(country_dirs), "country WKTs to use for GBIF download."))

# --- Selection for run ---
# Example: Process only the first country for a test run
country_dirs_to_process <- country_dirs[1] # Change index or remove [1] for full run

# --- Main Loop ---
for (country_dir in country_dirs_to_process) {

  country_name <- basename(country_dir)
  message(paste0("\n--- Processing Country: ", country_name, " ---"))

  # Note: This legacy script looks for a combined country WKT
  wkt_file_path <- file.path(country_dir, paste0("ramsar_sites_", country_name, ".wkt"))

  if (!file.exists(wkt_file_path)) {
    message(paste("  Skipping:", country_name, "- Combined WKT not found at:", wkt_file_path))
    next
  }

  wkt_string <- readLines(wkt_file_path, warn = FALSE)[1]

  # --- Initiate GBIF Download Request ---
  final_predicate <- get_gbif_predicates(wkt_string, gbif_base_filters)

  message(paste("  Requesting GBIF data for", country_name, "..."))
  download_key <- tryCatch({
    occ_download(
      final_predicate,
      format = "SIMPLE_CSV",
      user = getOption("gbif_user"),
      pwd = getOption("gbif_pwd"),
      email = getOption("gbif_email")
    )
  }, error = function(e) {
    message(paste("  Error requesting GBIF download for", country_name, ":", e$message))
    return(NULL)
  })

  if (is.null(download_key)) next

  message(paste("  Download initiated. Key:", download_key))

  # --- Wait for Download Completion and Retrieve Data ---
  message("  Waiting for download to complete...")
  download_status <- occ_download_wait(download_key)

  if (download_status$status != "SUCCEEDED") {
    message(paste("  Download failed. Status:", download_status$status))
    next
  }

  message(paste("  Download succeeded. Retrieving data..."))
  country_output_dir <- file.path(output_gbif_base_dir, country_name)
  if (!dir.exists(country_output_dir)) dir.create(country_output_dir, recursive = TRUE)

  zip_file_path <- occ_download_get(download_key, path = country_output_dir, overwrite = TRUE)
  
  # --- Unzip and Rename ---
  message(paste("  Unzipping..."))
  extracted_files <- utils::unzip(zip_file_path, exdir = country_output_dir)
  
  data_files <- extracted_files[grepl("\\.(csv|txt)$", extracted_files, ignore.case = TRUE)]

  if (length(data_files) > 0) {
    new_csv_path <- file.path(country_output_dir, paste0("gbif_occurrences_", country_name, ".csv"))
    file.rename(data_files[1], new_csv_path)
    message(paste("  Saved to:", new_csv_path))
  }

  if (file.exists(zip_file_path)) file.remove(zip_file_path)
}

message("\n--- GBIF data download process complete. ---")