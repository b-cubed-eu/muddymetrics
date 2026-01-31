# --- Setup ---
lapply(list.files("R", full.names = TRUE), source)
library(rgbif)
library(sf)
library(dplyr)
library(utils)

# --- Pre-flight Check ---
check_gbif_credentials()

# --- Configuration ---
wkt_input_base_dir <- "inst/extdata/ramsar_sites_wkt"
if (!dir.exists(wkt_input_base_dir)) {
  stop("WKT input directory not found: ", wkt_input_base_dir)
}

output_gbif_base_dir <- "inst/extdata/gbif_occurrences_by_ramsar_site"
if (!dir.exists(output_gbif_base_dir)) {
  dir.create(output_gbif_base_dir, recursive = TRUE)
}

gbif_base_filters <- list(
  hasGeospatialIssue = FALSE,
  hasCoordinate = TRUE,
  basisOfRecord = c("OBSERVATION", "HUMAN_OBSERVATION", "PRESERVED_SPECIMEN", "LITERATURE", "LIVING_SPECIMEN", "MATERIAL_SAMPLE", "MACHINE_OBSERVATION")
)

# --- Find all individual WKT files ---
all_wkt_files <- list.files(wkt_input_base_dir, pattern = "\\.wkt$", recursive = TRUE, full.names = TRUE)

if (length(all_wkt_files) == 0) {
  stop("No WKT files found in: ", wkt_input_base_dir)
}

message(paste("Found", length(all_wkt_files), "individual Ramsar site WKT files to process."))

# --- SELECT ONLY ONE SITE FOR TESTING ---
# Example: Process only the first site
wkt_files_to_process <- all_wkt_files[1] 

message(paste("\nRunning for:", basename(wkt_files_to_process)))

# --- Loop ---
for (wkt_filepath in wkt_files_to_process) {

  path_parts <- strsplit(wkt_filepath, .Platform$file.sep, fixed = TRUE)[[1]]
  country_name <- path_parts[length(path_parts) - 1]
  filename_full <- basename(wkt_filepath)
  sitename <- tools::file_path_sans_ext(filename_full)

  message(paste0("\n--- Processing Site: ", sitename, " (", country_name, ") ---"))

  wkt_string <- readLines(wkt_filepath, warn = FALSE)[1]

  # --- Initiate GBIF Download Request ---
  final_predicate <- get_gbif_predicates(wkt_string, gbif_base_filters)

  message(paste("  Requesting GBIF data..."))
  download_key <- tryCatch({
    occ_download(
      final_predicate,
      format = "SIMPLE_CSV",
      user = getOption("gbif_user"),
      pwd = getOption("gbif_pwd"),
      email = getOption("gbif_email")
    )
  }, error = function(e) {
    message(paste("  Error requesting download:", e$message))
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
  site_gbif_output_dir <- file.path(output_gbif_base_dir, country_name, sitename)
  if (!dir.exists(site_gbif_output_dir)) dir.create(site_gbif_output_dir, recursive = TRUE)

  zip_file_path <- occ_download_get(download_key, path = site_gbif_output_dir, overwrite = TRUE)
  
  # --- Unzip and Rename ---
  message(paste("  Unzipping..."))
  extracted_files <- utils::unzip(zip_file_path, exdir = site_gbif_output_dir)
  data_files <- extracted_files[grepl("\\.(csv|txt)$
", extracted_files, ignore.case = TRUE)]

  if (length(data_files) > 0) {
    new_csv_path <- file.path(site_gbif_output_dir, paste0("gbif_occurrences_", sitename, ".csv"))
    file.rename(data_files[1], new_csv_path)
    message(paste("  Saved to:", new_csv_path))
  }

  if (file.exists(zip_file_path)) file.remove(zip_file_path)
}

message("\n--- GBIF site-by-site download process complete. ---")