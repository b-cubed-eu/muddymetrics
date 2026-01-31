# --- Setup ---
lapply(list.files("R", full.names = TRUE), source)
library(rgbif)
library(sf)
library(dplyr)

# --- Pre-flight Check ---
check_gbif_credentials()

# --- Configuration ---
# Define the base directory where your BBOX WKT files are located
wkt_input_base_dir <- "inst/extdata/ramsar_sites_wkt"

# Directory to save the downloaded GBIF Data Cubes (CSV output)
output_gbif_cubes_base_dir <- "inst/extdata/gbif_data_cubes_mgrs_csv" 
if (!dir.exists(output_gbif_cubes_base_dir)) {
  dir.create(output_gbif_cubes_base_dir, recursive = TRUE)
}

# --- Cube specific parameters ---
mgrs_resolution_meters <- 1000 # 1km resolution
cube_measure <- "COUNT(DISTINCT speciesKey) AS unique_species_count"

# --- SELECT A SINGLE RAMSAR SITE WKT FILE FOR TESTING ---
all_wkt_files <- list.files(wkt_input_base_dir, pattern = "\.wkt$", recursive = TRUE, full.names = TRUE)
if (length(all_wkt_files) == 0) stop("No WKT files found.")

test_wkt_filepath <- all_wkt_files[1]

message(paste0("\n--- Requesting GBIF Data Cube (MGRS ", mgrs_resolution_meters/1000, "km, CSV) for: ", basename(test_wkt_filepath), " ---"))

# 1. Read the WKT content and extract coordinates
wkt_string <- readLines(test_wkt_filepath, warn = FALSE)[1]
bbox_sf <- st_as_sfc(wkt_string, crs = 4326)
bbox_coords <- st_bbox(bbox_sf)

min_lat <- bbox_coords["ymin"]
max_lat <- bbox_coords["ymax"]
min_lon <- bbox_coords["xmin"]
max_lon <- bbox_coords["xmax"]

# 2. Construct the SQL query
sql_query <- sprintf(
  "\n  SELECT\n    GBIF_MGRSCode(%d, decimalLatitude, decimalLongitude, 0.) AS mgrs_grid_cell,\n    %s\n  FROM\n    occurrence\n  WHERE\n    decimalLatitude >= %f AND decimalLatitude <= %f\n    AND decimalLongitude >= %f AND decimalLongitude <= %f\n    AND hasCoordinate = TRUE\n    AND coordinateUncertaintyInMeters <= 50000\n    AND countryCode IS NOT NULL\n  GROUP BY\n    mgrs_grid_cell\n",
  mgrs_resolution_meters, cube_measure, min_lat, max_lat, min_lon, max_lon
)

message("\nGenerated SQL Query:")
message(sql_query)

# 3. Initiate the download
message("Initiating GBIF Data Cube (CSV) download via SQL API...")
download_object <- tryCatch({
  occ_download_sql(
    q = sql_query,
    user = getOption("gbif_user"),
    pwd = getOption("gbif_pwd"),
    email = getOption("gbif_email")
  )
}, error = function(e) {
  message(paste("Error initiating SQL download:", e$message))
  return(NULL)
})

if (is.null(download_object)) stop("Failed to initiate download.")

download_key <- download_object[1]
message(paste("Request submitted. Key:", download_key))

# 4. Wait for the download to complete
message("Waiting for data cube to be processed...")
gbif_download_status <- occ_download_wait(x = download_key)

if (gbif_download_status$status != "SUCCEEDED") {
  stop(paste("Download failed. Status:", gbif_download_status$status))
}

# 5. Retrieve and Organize
path_parts <- strsplit(test_wkt_filepath, .Platform$file.sep, fixed = TRUE)[[1]]
country_name <- path_parts[length(path_parts) - 1]
sitename <- tools::file_path_sans_ext(basename(test_wkt_filepath))

site_cube_output_dir <- file.path(output_gbif_cubes_base_dir, country_name, sitename)
if (!dir.exists(site_cube_output_dir)) dir.create(site_cube_output_dir, recursive = TRUE)

output_csv_path <- file.path(site_cube_output_dir, paste0("gbif_datacube_", sitename, ".csv"))

message("Retrieving data...")
downloaded_zip_path <- occ_download_get(key = download_key, path = site_cube_output_dir, overwrite = TRUE)

# Unzip and rename
extracted_files <- utils::unzip(downloaded_zip_path, exdir = site_cube_output_dir)
potential_data_files <- grep("\.(txt|csv)$", extracted_files, ignore.case = TRUE, value = TRUE)

if (length(potential_data_files) > 0) {
  file_sizes <- file.info(potential_data_files)$size
  csv_file <- potential_data_files[which.max(file_sizes)]
  file.rename(csv_file, output_csv_path)
  message(paste("Saved to:", output_csv_path))
}

if (file.exists(downloaded_zip_path)) file.remove(downloaded_zip_path)

message("\n--- GBIF Data Cube download complete. ---")