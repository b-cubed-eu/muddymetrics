library(sf)
library(dplyr)
library(rgbif) # Assuming you'll use this later for GBIF downloads
library(utils) # For unzip if you combine with GBIF download part

# --- 1. Load and prepare global Ramsar shapefile ---
message("Loading global Ramsar shapefile...")
ramsar_global <- sf::st_read("inst/extdata/Ramsar_boundaries/features_publishedPolygon.shp")

ramsar_global %>% # fix bad country names
  dplyr::mutate(
    country_en = dplyr::case_when(
      country_en == "C\xf4te d'Ivoire" ~ "Cote dIvoire",
      country_en == "T\xfcrkiye" ~ "Turkey",
      country_en == "Viet Nam" ~ "Vietnam",
      TRUE ~ country_en
    ),
    ramsar_site_id = paste0("site_", row_number()) # Creates unique ID like "site_1", "site_2"
  ) -> ramsar_global

message(paste("Loaded", nrow(ramsar_global), "individual Ramsar sites, and created unique IDs."))

# --- 2. Configuration for output directories ---
output_wkt_base_dir <- "inst/extdata/ramsar_sites_wkt"
if (!dir.exists(output_wkt_base_dir)) {
  dir.create(output_wkt_base_dir, recursive = TRUE)
  message(paste("Created base WKT output directory:", output_wkt_base_dir))
}

# --- 3. Process each individual Ramsar site ---
message("\nStarting WKT conversion for individual Ramsar sites...")
for (j in 1:nrow(ramsar_global)) {
  ramsar_site <- ramsar_global[j, ]

  ramsar_site_unique_id <- ramsar_site$ramsar_site_id
  site_name_original <- ramsar_site$officialna
  country_name <- ramsar_site$country_en

  # --- Robust Sanitization of site name for file paths ---
  site_name_base <- if (is.na(site_name_original) || length(site_name_original) == 0 || site_name_original == "") {
    as.character(ramsar_site_unique_id)
  } else {
    as.character(site_name_original)
  }

  site_name_base_utf8 <- iconv(site_name_base, from = "", to = "UTF-8", sub = "_")

  site_name_sanitized <- gsub("[[:punct:]\\s]+", "_", site_name_base_utf8)
  site_name_sanitized <- gsub("_+", "_", site_name_sanitized)
  site_name_sanitized <- gsub("^_|_$", "", site_name_sanitized)

  if (nchar(site_name_sanitized) == 0) {
    site_name_sanitized <- as.character(ramsar_site_unique_id)
  }
  # --- End Robust Sanitization ---

  message(paste0("\n--- Processing site: ", site_name_original, " (ID: ", ramsar_site_unique_id, ", Country: ", country_name, ") ---"))

  # Create country-specific subdirectory
  country_wkt_dir <- file.path(output_wkt_base_dir, country_name)
  if (!dir.exists(country_wkt_dir)) {
    dir.create(country_wkt_dir, recursive = TRUE)
    message(paste("  Created country directory:", country_wkt_dir))
  }

  # Define output file path for WKT for this specific site
  output_wkt_path <- file.path(country_wkt_dir, paste0(ramsar_site_unique_id, "_", site_name_sanitized, ".wkt"))

  # --- Geometry Operations ---
  # Remove the explicit sf_use_s2(FALSE) call here.
  # This means st_make_valid will run with the default s2 setting (which is TRUE for sf >= 1.0.0).
  # We still capture and restore the original setting, though it's less critical now.
  original_s2_setting <- sf_use_s2() # Capture original setting

  ramsar_site_valid <- tryCatch({
    st_make_valid(ramsar_site)
  }, error = function(e) {
    message(paste("  Error making geometry valid for site", ramsar_site_unique_id, ":", e$message))
    return(NULL)
  }, finally = {
    sf_use_s2(original_s2_setting) # Always restore s2 setting
  })

  if (is.null(ramsar_site_valid) || st_is_empty(ramsar_site_valid)) {
    message(paste("  Skipping WKT creation for site", ramsar_site_unique_id, "due to invalid or empty geometry."))
    next
  }

  # --- Convert the site geometry to WKT ---
  wkt_string <- tryCatch({
    st_as_text(ramsar_site_valid$geom) # Assuming 'geom' is your geometry column
  }, error = function(e) {
    message(paste("  Error converting site", ramsar_site_unique_id, "to WKT:", e$message))
    return(NULL)
  })

  if (is.null(wkt_string) || nchar(wkt_string) < 10) {
    message(paste("  Skipping site", ramsar_site_unique_id, "- WKT string is empty or invalid."))
    next
  }

  # --- Save the WKT string to a .wkt file ---
  message(paste("  Saving WKT to:", output_wkt_path))
  writeLines(wkt_string, output_wkt_path)
}

message("\n--- Individual Ramsar site WKT conversion complete. ---")




message("\n--- Quick WKT File Verification ---")

# 1. Get the total number of Ramsar sites that were processed (or attempted)
total_expected_wkt_files <- nrow(ramsar_global)

# 2. Count the number of .wkt files actually created on disk
#    list.files will search recursively in the output_wkt_base_dir
#    pattern = "\\.wkt$" ensures we only count files ending with .wkt
actual_wkt_files_on_disk <- length(list.files(output_wkt_base_dir, pattern = "\\.wkt$", recursive = TRUE))

message(paste("Expected WKT files (based on input sites):", total_expected_wkt_files))
message(paste("Actual WKT files found on disk:", actual_wkt_files_on_disk))

if (actual_wkt_files_on_disk == total_expected_wkt_files) {
  message("SUCCESS: All expected WKT files appear to have been created.")
} else if (actual_wkt_files_on_disk < total_expected_wkt_files) {
  message(paste("WARNING: Only", actual_wkt_files_on_disk, "out of", total_expected_wkt_files, "WKT files were found."))
  message("This indicates some sites might have failed or been skipped.")
  # To find out which ones, you'd need the more detailed logging previously provided.
} else {
  message(paste("INFO: More WKT files found on disk (", actual_wkt_files_on_disk, ") than expected (", total_expected_wkt_files, ")."))
  message("This might be due to previous runs leaving files behind.")
}




message("\n--- Detailed WKT File Verification ---")

# 1. Get the total number of Ramsar sites that were processed (or attempted)
total_expected_wkt_files <- nrow(ramsar_global)

# 2. Get the list of all unique ramsar_site_id values that should have a WKT file
expected_ramsar_ids <- ramsar_global$ramsar_site_id

# 3. Count the number of .wkt files actually created on disk
generated_wkt_files_paths <- list.files(output_wkt_base_dir, pattern = "\\.wkt$", recursive = TRUE, full.names = TRUE)
actual_wkt_files_on_disk_count <- length(generated_wkt_files_paths)

message(paste("Expected WKT files (based on input sites):", total_expected_wkt_files))
message(paste("Actual WKT files found on disk:", actual_wkt_files_on_disk_count))

if (actual_wkt_files_on_disk_count == total_expected_wkt_files) {
  message("SUCCESS: All expected WKT files appear to have been created.")
} else if (actual_wkt_files_on_disk_count < total_expected_wkt_files) {
  message(paste("WARNING: Only", actual_wkt_files_on_disk_count, "out of", total_expected_wkt_files, "WKT files were found."))
  message("Identifying the missing file(s)...")

  # Extract ramsar_site_id from the filenames actually found on disk
  # This relies on the naming convention: "ramsar_site_id_sanitized_site_name.wkt"
  # We take the first part split by underscore.
  actual_ramsar_ids_on_disk <- sapply(basename(generated_wkt_files_paths), function(filename) {
    strsplit(filename, "_", fixed = TRUE)[[1]][1]
  })
  actual_ramsar_ids_on_disk <- unique(actual_ramsar_ids_on_disk) # Ensure uniqueness in case of anomalies

  # Find which IDs are in the 'expected' list but not in the 'actual' list
  missing_ramsar_ids <- setdiff(expected_ramsar_ids, actual_ramsar_ids_on_disk)

  if (length(missing_ramsar_ids) > 0) {
    message("\n--- Details of Missing WKT File(s) ---")
    missing_site_details <- ramsar_global %>%
      filter(ramsar_site_id %in% missing_ramsar_ids) %>%
      select(ramsar_site_id, officialna, country_en) # Select relevant columns for identification

    if (nrow(missing_site_details) > 0) {
      print(missing_site_details)
    } else {
      message("Could not find details for the missing IDs in the original data. This is unexpected.")
    }
    message("\nThese sites likely failed during WKT conversion or file writing.")

  } else {
    message("No specific missing IDs found despite count mismatch. This is unusual and might indicate a naming or parsing issue.")
  }

} else { # actual_wkt_files_on_disk_count > total_expected_wkt_files
  message(paste("INFO: More WKT files found on disk (", actual_wkt_files_on_disk_count, ") than expected (", total_expected_wkt_files, ")."))
  message("This might be due to previous runs leaving files behind, or a discrepancy in file counting/naming.")
  # To identify these "extra" files, you could reverse the setdiff:
  # extra_files_on_disk <- setdiff(actual_ramsar_ids_on_disk, expected_ramsar_ids)
  # message("Extra files found (IDs):", paste(extra_files_on_disk, collapse = ", "))
}




# --- 3. Initialize a data frame to store processing status ---
site_processing_status <- data.frame(
  ramsar_site_id = character(),
  officialna = character(),
  country_en = character(),
  status = character(),
  message = character(),
  stringsAsFactors = FALSE
)

# --- 4. SELECT THE SINGLE RAMSAR SITE TO TEST ---
# IMPORTANT: Use the ID of the site that is consistently failing with GBIF.
single_test_site_id <- "site_1" # <--- REPLACE WITH THE ID OF THE SITE YOU WANT TO TEST!

message(paste("\n*** Processing ONLY the single test site with ID:", single_test_site_id, "***"))

# --- Loop, but only process the selected site ---
site_row_index <- which(ramsar_global$ramsar_site_id == single_test_site_id)

if (length(site_row_index) == 0) {
  stop(paste("Site with ID", single_test_site_id, "not found in the loaded Ramsar data. Please check the ID."))
}

for (j in site_row_index) { # Loop will run only for the specified row index

  ramsar_site <- ramsar_global[j, ]

  ramsar_site_unique_id <- ramsar_site$ramsar_site_id
  site_name_original <- ramsar_site$officialna
  country_name <- ramsar_site$country_en

  # --- Start a comprehensive tryCatch block for each site's processing ---
  processing_result <- tryCatch({
    # --- Robust Sanitization of site name for file paths ---
    site_name_base <- if (is.na(site_name_original) || length(site_name_original) == 0 || site_name_original == "") {
      as.character(ramsar_site_unique_id)
    } else {
      as.character(site_name_original)
    }

    site_name_base_utf8 <- iconv(site_name_base, from = "", to = "UTF-8", sub = "_")
    if (is.na(site_name_base_utf8)) {
      stop("Encoding conversion failed for site name.")
    }

    site_name_sanitized <- gsub("[[:punct:]\\s]+", "_", site_name_base_utf8)
    site_name_sanitized <- gsub("_+", "_", site_name_sanitized)
    site_name_sanitized <- gsub("^_|_$", "", site_name_sanitized)
    if (nchar(site_name_sanitized) == 0) {
      site_name_sanitized <- as.character(ramsar_site_unique_id)
    }
    # --- End Robust Sanitization ---

    message(paste0("\n--- Processing site: ", site_name_original, " (ID: ", ramsar_site_unique_id, ", Country: ", country_name, ") ---"))
    message("  (NOTE: For this test, we are generating a BOUNDING BOX WKT, not the exact shape.)")
    message("  This will cover a potentially larger area for GBIF download.")

    # Create country-specific subdirectory
    country_wkt_dir <- file.path(output_wkt_base_dir, country_name)
    if (!dir.exists(country_wkt_dir)) {
      dir.create(country_wkt_dir, recursive = TRUE)
      message(paste("  Created country directory:", country_wkt_dir))
    }

    # IMPORTANT: Adjust filename for clarity that it's a bounding box
    output_wkt_path <- file.path(country_wkt_dir, paste0(ramsar_site_unique_id, "_", site_name_sanitized, "_bbox.wkt")) # <--- Added "_bbox"

    # --- Geometry Operations: Create Bounding Box ---
    # Get the bounding box coordinates
    bbox <- st_bbox(ramsar_site)

    # Convert the bounding box to an sf polygon object (a rectangle)
    # This automatically creates a valid POLYGON WKT
    bbox_polygon <- st_as_sfc(bbox)

    # Ensure it has the correct CRS (WGS84, EPSG:4326)
    if (st_crs(bbox_polygon) != st_crs(4326)) {
      bbox_polygon <- st_transform(bbox_polygon, 4326)
    }

    # Convert the bounding box polygon to WKT
    wkt_string <- st_as_text(bbox_polygon)

    if (is.null(wkt_string) || nchar(wkt_string) < 10) {
      stop("WKT string is empty or too short after bounding box conversion.")
    }

    # Save the WKT string to a .wkt file
    writeLines(wkt_string, output_wkt_path)

    list(status = "SUCCESS", message = "Bounding box WKT created successfully.", wkt_path = output_wkt_path)

  }, error = function(e) {
    list(status = "FAILED", message = paste("Error:", e$message))
  }, finally = {
    # No s2 setting restore needed as we are directly creating a new sf object from bbox
    # and not modifying existing complex geometry that might have s2 issues.
  })

  site_processing_status <- rbind(site_processing_status, data.frame(
    ramsar_site_id = ramsar_site_unique_id,
    officialna = ifelse(is.null(site_name_original) || is.na(site_name_original), NA_character_, site_name_original),
    country_en = country_name,
    status = processing_result$status,
    message = processing_result$message,
    stringsAsFactors = FALSE
  ))

  if (processing_result$status == "SUCCESS") {
    message(paste("  WKT saved to:", processing_result$wkt_path))
  } else {
    message(paste("  Processing Failed/Skipped for site", ramsar_site_unique_id, ":", processing_result$message))
  }
}

message("\n--- Single Ramsar site BOUNDING BOX WKT conversion test complete. ---")

# --- Report for the single site processed ---
message("\n--- Test Site Processing Result ---")
print(site_processing_status)
