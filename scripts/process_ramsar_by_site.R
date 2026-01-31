# This script processes the global Ramsar sites shapefile, consolidates geometries, and generates WKT files for each site.
# It also includes robust error handling and sanitization for site names.

# --- 1. Install and load required packages ---
packages <- c("sf", "dplyr", "stringi") # Add all packages your script uses
for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}

# --- 2. Configure paths ---
ramsar_global_path <- "inst/extdata/Ramsar_boundaries/features_publishedPolygon.shp"
output_wkt_base_dir <- "inst/extdata/ramsar_sites_wkt"

# --- 3. Load and prepare global Ramsar shapefile ---
message("Loading global Ramsar shapefile...")
ramsar_global <- sf::st_read(ramsar_global_path, quiet = TRUE)
if (inherits(ramsar_global, "try-error")) {
  stop("FATAL ERROR: Could not read the global Ramsar shapefile. Please check the file path and format.")
}
ramsar_global %>% # fix bad country names
  dplyr::mutate(
    country_en = dplyr::case_when(
      country_en == "C\xf4te d'Ivoire" ~ "Cote dIvoire",
      country_en == "T\xfcrkiye" ~ "Turkey",
      country_en == "Viet Nam" ~ "Vietnam",
      TRUE ~ country_en
    ),
    ramsar_site_id = paste0("site_", ramsarid) # Creates unique ID like "site_1", "site_2"
  ) -> ramsar_global

message(paste("Loaded", nrow(ramsar_global), "individual Ramsar sites, and created unique IDs."))

# --- 4. Create output directory if missing ---
if (!dir.exists(output_wkt_base_dir)) {
  dir.create(output_wkt_base_dir, recursive = TRUE)
  message(paste("Created base WKT output directory:", output_wkt_base_dir))
}

# --- 5. Consolidate duplicate ramsar sites ---
# Check for duplicate ramsar_site_id values BEFORE consolidation
duplicate_ids_before_consolidation <- ramsar_global$ramsar_site_id[duplicated(ramsar_global$ramsar_site_id)]
if (length(duplicate_ids_before_consolidation) > 0) {
  message("\n*** WARNING: Duplicate 'ramsar_site_id' values found BEFORE consolidation! ***")
  message("This means some sites are represented by multiple features in the original shapefile.")
  message("Duplicate IDs (first 10):", paste(head(unique(duplicate_ids_before_consolidation), 10), collapse = ", "))
} else {
  message("\n--- No duplicate 'ramsar_site_id' values found BEFORE consolidation. ---")
}

# Consolidate geometries for sites with the same ramsar_site_id
message("\n--- Consolidating geometries for duplicate Ramsar Site IDs... ---")

# Apply st_make_valid to the entire ramsar_global object BEFORE grouping and unioning
message("  [Timing] Starting st_make_valid on entire ramsar_global object:", Sys.time())
ramsar_global_valid <- tryCatch({
  st_make_valid(ramsar_global)
}, error = function(e) {
  stop("FATAL ERROR: Could not make all global Ramsar geometries valid before consolidation. Error: ", e$message)
})
message("  [Timing] Finished st_make_valid on entire ramsar_global object:", Sys.time())


# Ensure the sf_use_s2 setting is captured and restored
original_s2_setting <- sf_use_s2() # Capture original setting
sf_use_s2(FALSE) # Disable s2 for this operation to avoid issues with complex geometries

# Now, group and union the *valid* geometries with a robust error handling for st_union
ramsar_global_consolidated <- ramsar_global_valid %>%
  group_by(ramsar_site_id, country_en, officialna) %>% # Include officialna to retain it for sanitization
  summarise(
    geometry = {
      # Attempt to union the geometries for the current group
      union_result <- tryCatch({
        st_union(geometry)
      }, error = function(e) {
        # If st_union fails, log the error and return a GEOMETRYCOLLECTION of the individual parts
        message(paste0("    [Consolidation Error] Failed to union geometry for site ID: ", first(ramsar_site_id), " (Official Name: ", first(officialna), "). Error: ", e$message, ". Returning individual valid geometries as GEOMETRYCOLLECTION."))
        # st_combine creates a single multi-geometry object without topological checks
        # st_cast ensures it's explicitly a GEOMETRYCOLLECTION
        st_cast(st_combine(geometry), "GEOMETRYCOLLECTION")
      })
      union_result
    },
    .groups = "drop"
  ) %>%
  st_as_sf()

sf_use_s2(original_s2_setting) # Always restore s2 setting

message(paste("Consolidated from", nrow(ramsar_global), "features to", nrow(ramsar_global_consolidated), "unique Ramsar sites."))
message("\n--- Details of first 5 CONSOLIDATED Ramsar sites (if available): ---")
print(head(ramsar_global_consolidated %>% select(ramsar_site_id, officialna, country_en), 5))

# --- 5. Process each individual Ramsar site ---
message("\nStarting WKT conversion for individual Ramsar sites...")
for (j in 1:nrow(ramsar_global_consolidated)) {
  ramsar_site <- ramsar_global_consolidated[j, ]

  ramsar_site_unique_id <- ramsar_site$ramsar_site_id
  site_name_original <- ramsar_site$officialna
  country_name <- ramsar_site$country_en

  # --- Robust Sanitization of site name for file paths ---
  site_name_base <- if (is.na(site_name_original) || length(site_name_original) == 0 || site_name_original == "") {
    as.character(ramsar_site_unique_id)
  } else {
    as.character(site_name_original)
  }

  site_name_base_utf8 <- iconv(site_name_base, from = "Latin1", to = "UTF-8", sub = "?")

  site_name_transliterated <- stringi::stri_trans_general(site_name_base_utf8, "Latin-ASCII")

  site_name_no_apostrophe <- gsub("'", "", site_name_transliterated, fixed = TRUE)

  site_name_sanitized <- gsub("[^[:alnum:]\\.\\_\\-]+", "_", site_name_no_apostrophe)
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
  ramsar_site_valid <- tryCatch({
    st_make_valid(ramsar_site)
  }, error = function(e) {
    message(paste("  Error making geometry valid for site", ramsar_site_unique_id, ":", e$message))
    return(NULL)
  })

  if (is.null(ramsar_site_valid) || st_is_empty(ramsar_site_valid)) {
    message(paste("  Skipping WKT creation for site", ramsar_site_unique_id, "due to invalid or empty geometry."))
    next
  }

  # --- Convert the site geometry to WKT ---
  wkt_string <- tryCatch({
    st_as_text(ramsar_site_valid$geometry) # Assuming 'geom' is your geometry column
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


# --- 6. Verification of WKT Files ---
message("\n--- Detailed WKT File Verification ---")

# Get the total number of Ramsar sites that were processed (or attempted)
total_expected_wkt_files <- nrow(ramsar_global_consolidated)

# Get the list of all unique ramsar_site_id values that should have a WKT file
expected_ramsar_ids <- ramsar_global_consolidated$ramsar_site_id

# Count the number of .wkt files actually created on disk
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

