library(rgbif)
library(sf)
library(dplyr)

# --- Configuration for GBIF downloads ---
# Make sure your GBIF credentials are set as R options or environment variables!
options(gbif_user = "shawn.dove", gbif_pwd = "Iliketoast43@!S", gbif_email = "shawn.dove@hotmail.com")
# Replace with your actual credentials
gbif_user <- getOption("gbif_user")
gbif_pwd <- getOption("gbif_pwd")
gbif_email <- getOption("gbif_email")

if (is.null(gbif_user) || is.null(gbif_pwd) || is.null(gbif_email)) {
  stop("GBIF credentials (user, pwd, email) are not set. Please set them using options() or environment variables before running.")
}

# Define the base directory where your BBOX WKT files are located
wkt_input_base_dir <- "inst/extdata/ramsar_sites_wkt_bbox"

# Directory to save the downloaded GBIF Data Cubes (CSV output)
output_gbif_cubes_base_dir <- "inst/extdata/gbif_data_cubes_mgrs_csv" # Output will be CSV
if (!dir.exists(output_gbif_cubes_base_dir)) {
  dir.create(output_gbif_cubes_base_dir, recursive = TRUE)
  message(paste("Created base GBIF data cube CSV output directory:", output_gbif_cubes_base_dir))
}

# --- Cube specific parameters ---
mgrs_resolution_meters <- 1000 # 1km resolution (use 10000 for 10km, 100 for 100m etc.)
cube_measure <- "COUNT(DISTINCT speciesKey) AS unique_species_count" # Or COUNT(*) AS occurrence_count

# --- SELECT A SINGLE RAMSAR SITE BBOX WKT FILE FOR TESTING ---
# IMPORTANT: Adjust this path to one of your successfully generated _bbox.wkt files.
# DOUBLE CHECK THIS PATH AND THE FILE'S EXISTENCE!
test_bbox_wkt_filepath <- "inst/extdata/ramsar_sites_wkt/Switzerland/site_1_Le Rh_ne genevoi_ _ Vallon_ de l_Allondon et de la Laire_bbox.wkt"

if (!file.exists(test_bbox_wkt_filepath)) {
  stop("The specified test_bbox_wkt_filepath does not exist: ", test_bbox_wkt_filepath,
       ". Please choose a valid path to an existing BBOX WKT file.")
}

message(paste0("\n--- Requesting GBIF Data Cube (MGRS ", mgrs_resolution_meters/1000, "km, CSV) for: ", basename(test_bbox_wkt_filepath), " ---"))

# 1. Read the BBOX WKT content as plain text, then convert to sf object to extract coordinates
message(paste("Reading BBOX WKT content from:", test_bbox_wkt_filepath, "to extract coordinates."))
wkt_string_from_file <- readLines(test_bbox_wkt_filepath, warn = FALSE)[1]

if (is.null(wkt_string_from_file) || nchar(wkt_string_from_file) < 10) {
  stop("Failed to read a valid WKT string from: ", test_bbox_wkt_filepath)
}

bbox_sf <- st_as_sfc(wkt_string_from_file, crs = 4326) # Explicitly set CRS to WGS84
if (st_crs(bbox_sf) != st_crs(4326)) {
  bbox_sf <- st_transform(bbox_sf, 4326) # Ensure it's WGS84
}
bbox_coords <- st_bbox(bbox_sf)

min_lat <- bbox_coords["ymin"]
max_lat <- bbox_coords["ymax"]
min_lon <- bbox_coords["xmin"]
max_lon <- bbox_coords["xmax"]

message(paste0("  Extracted bounds: Lat (", min_lat, " to ", max_lat, "), Lon (", min_lon, " to ", max_lon, ")"))


# 2. Construct the SQL query with GBIF_MGRSCode and simple bounding box filter
# Removed the HAVING clause as it's not supported by the GBIF SQL API.
sql_query <- sprintf("
  SELECT
    GBIF_MGRSCode(%d, decimalLatitude, decimalLongitude, 0.) AS mgrs_grid_cell,
    %s
  FROM
    occurrence
  WHERE
    decimalLatitude >= %f AND decimalLatitude <= %f
    AND decimalLongitude >= %f AND decimalLongitude <= %f
    AND hasCoordinate = TRUE
    AND coordinateUncertaintyInMeters <= 50000 -- Optional: filter for more precise coords
    AND countryCode IS NOT NULL -- Essential for most analyses
  GROUP BY
    mgrs_grid_cell
", mgrs_resolution_meters, cube_measure, min_lat, max_lat, min_lon, max_lon) # Inserting the extracted coordinates

message("\nGenerated SQL Query for GBIF CSV Cube:")
message(sql_query)

# 3. Initiate the download using occ_download_sql() with 'q'
message("Initiating GBIF Data Cube (CSV) download via SQL API...")
download_object <- tryCatch({
  occ_download_sql(
    q = sql_query, # Changed to 'q' based on your feedback
    user = gbif_user,
    pwd = gbif_pwd,
    email = gbif_email
  )
}, error = function(e) {
  message(paste("Error initiating SQL download:", e$message))
  return(NULL)
})

if (is.null(download_object)) {
  stop("Failed to initiate download. Check SQL query and GBIF credentials.")
}

download_key <- download_object[1]
message(paste("GBIF Data Cube (CSV) request submitted. Download Key:", download_key))
message("Waiting for data cube to be processed by GBIF... (This may take significant time)")

# 4. Wait for the download to complete
gbif_download_status <- tryCatch({
  occ_download_wait(x = download_key)
}, error = function(e) {
  message(paste("Error waiting for GBIF Data Cube download:", e$message))
  return(NULL)
})

if (is.null(gbif_download_status) || gbif_download_status$status != "SUCCEEDED") {
  stop(paste("GBIF Data Cube (CSV) download did not succeed for key:", download_key, ". Status:", gbif_download_status$status,
             ". Check GBIF portal for more details on download", download_key))
}

message(paste("GBIF Data Cube (CSV) download for key", download_key, "SUCCEEDED."))

# 5. Extract site ID and name for output organization
filename_parts <- strsplit(basename(test_bbox_wkt_filepath), "_", fixed = TRUE)[[1]]
ramsar_site_id_test <- filename_parts[1]
sanitized_site_name_test <- paste(filename_parts[2:(length(filename_parts)-2)], collapse = "_")
test_country_name <- basename(dirname(test_bbox_wkt_filepath))

# Create specific output directory for this site's GBIF data cube
site_cube_output_dir <- file.path(output_gbif_cubes_base_dir, test_country_name, paste0(ramsar_site_id_test, "_", sanitized_site_name_test))
if (!dir.exists(site_cube_output_dir)) {
  dir.create(site_cube_output_dir, recursive = TRUE)
  message(paste("Created GBIF data cube CSV output directory for test site:", site_cube_output_dir))
}

# Define the output CSV filename
output_csv_filename <- paste0("gbif_datacube_", ramsar_site_id_test, "_MGRS_", mgrs_resolution_meters/1000, "km.csv")
output_csv_path <- file.path(site_cube_output_dir, output_csv_filename)

message("Retrieving GBIF Data Cube (CSV)...")
downloaded_zip_path <- tryCatch({
  occ_download_get(key = download_key, path = site_cube_output_dir)
}, error = function(e) {
  message(paste("Error retrieving GBIF Data Cube (CSV):", e$message))
  return(NULL)
})

if (is.null(downloaded_zip_path) || !file.exists(downloaded_zip_path)) {
  stop("Failed to retrieve GBIF Data Cube (CSV). Check error messages above.")
}

message(paste("Downloaded zip file to:", downloaded_zip_path))

# Unzip and rename the main data file (SQL downloads are often 'occurrence.txt')
extracted_files <- utils::unzip(downloaded_zip_path, exdir = site_cube_output_dir)
# --- CORRECTED PART: Identify the main CSV/TXT data file ---
# Find all text or CSV files in the extracted list
potential_data_files <- grep("\\.(txt|csv)$", extracted_files, ignore.case = TRUE, value = TRUE)

if (length(potential_data_files) == 0) {
  stop("Could not find any .txt or .csv file in the downloaded zip. Extracted files: ", paste(basename(extracted_files), collapse = ", "))
} else if (length(potential_data_files) == 1) {
  # If only one potential data file, use it
  csv_file <- potential_data_files[1]
} else {
  # If multiple .txt/.csv files, assume the largest one is the main data file
  # This is common as GBIF zips often include metadata like 'citations.txt' or 'README.txt'
  file_sizes <- file.info(potential_data_files)$size
  csv_file <- potential_data_files[which.max(file_sizes)]
  message(paste("Multiple .txt/.csv files found. Selecting the largest one as the main data file:", basename(csv_file)))
}
# --- END CORRECTED PART ---

file.rename(csv_file[1], output_csv_path)
message(paste("GBIF Data Cube (CSV) saved to:", output_csv_path))

# (Optional) Clean up other extracted files and the original zip
other_files <- setdiff(extracted_files, csv_file[1])
if(length(other_files) > 0) file.remove(other_files)
file.remove(downloaded_zip_path)

message("\n--- Single site GBIF Pre-gridded Data Cube (MGRS CSV) download test complete. ---")
message(paste0("Your MGRS data cube in CSV format is in: ", output_csv_path))
message("This CSV should contain 'mgrs_grid_cell' and your chosen measure (e.g., 'unique_species_count').")
