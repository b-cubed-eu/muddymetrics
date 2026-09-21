# ============================================================================
# SCRIPT: Safe & Targeted Split GBIF MGRS Data Cubes for Remaining Continents (Resumed & Optimized)
# ============================================================================

library(sf)
library(data.table)
library(dplyr)
library(stringi)
library(mgrs)

# --- Configuration ---
old_shp_path <- "inst/extdata/Ramsar_boundaries/features_publishedPolygon.shp"
new_shp_path <- "inst/extdata/Ramsar_boundaries_160626/features_publishedPolygon.shp"
ramsar_wkt_base_dir <- "inst/extdata/ramsar_sites_wkt"
mgrs_column_name <- "mgrscellcode"
chunk_size <- 10000000 # 10 million rows per chunk (highly stable for Oceania/South America)
continents <- c("oceania", "southamerica")

# Resume configuration (where to start reading each continent)
resume_rows_by_cont <- list(
  oceania = 0,
  southamerica = 0
)

# --- Step 1: Identify target sites ---
cat("\n=== Identifying target sites ===\n")
old_shp <- st_read(old_shp_path, quiet = TRUE)
new_shp <- st_read(new_shp_path, quiet = TRUE)

old_ids <- sort(unique(old_shp$ramsarid))
new_ids <- sort(unique(new_shp$ramsarid))
common_ids <- intersect(old_ids, new_ids)

added_in_new <- setdiff(new_ids, old_ids)
geom_changed_sites <- c("2128", "2324", "514", "1110", "1384", "1386", "1681", "2126", "2127")

attr_changed_sites <- c()
for (rid in common_ids) {
  old_sub <- old_shp[old_shp$ramsarid == rid, ]
  new_sub <- new_shp[new_shp$ramsarid == rid, ]
  if (old_sub$officialna[1] != new_sub$officialna[1] ||
      old_sub$iso3[1] != new_sub$iso3[1] ||
      old_sub$country_en[1] != new_sub$country_en[1] ||
      old_sub$area_off[1] != new_sub$area_off[1]) {
    attr_changed_sites <- c(attr_changed_sites, rid)
  }
}

target_ramsarids <- unique(c(added_in_new, geom_changed_sites, attr_changed_sites))
cat("Total target Ramsar IDs to process: ", length(target_ramsarids), "\n")

# --- Step 2: Move existing CSV files of target sites to temp_holding (if not resuming) ---
cat("\n=== Moving existing CSV files of target sites to temp holding ===\n")
if (!dir.exists("temp_holding")) {
  dir.create("temp_holding", recursive = TRUE)
}

for (cont in continents) {
  # If we are resuming, skip clearing the existing partial files for this continent
  if (resume_rows_by_cont[[cont]] > 0) {
    cat(sprintf("Continent %s is resuming from line %s. Preserving existing partial CSV files.\n", 
                cont, format(resume_rows_by_cont[[cont]] + 2, scientific=FALSE)))
    next
  }
  
  data_dir <- file.path("inst/extdata", paste0("ramsar_site_data_100m_", cont))
  if (dir.exists(data_dir)) {
    for (fid in target_ramsarids) {
      csv_files <- list.files(data_dir, pattern = paste0("^site_", fid, "_.*_data\\.csv$"), recursive = TRUE, full.names = TRUE)
      for (csv_f in csv_files) {
        dest <- file.path("temp_holding", paste0(cont, "_", basename(csv_f)))
        # Avoid overwriting if it's already in temp_holding
        if (!file.exists(dest)) {
          success <- file.rename(csv_f, dest)
          if (success) {
            cat("Moved existing CSV:", csv_f, "->", dest, "\n")
          }
        } else {
          file.remove(csv_f)
          cat("Removed redundant CSV (already backed up in temp_holding):", csv_f, "\n")
        }
      }
    }
  }
}

# --- Step 3: Load target WKT files ---
cat("\n=== Loading target WKT files ===\n")

# Setup country to continent mapping
continents_list <- c("africa", "antarctica", "asia", "europe", "northamerica", "oceania", "southamerica")
country_to_cont <- list()
for (c_name in continents_list) {
  data_dir <- file.path("inst/extdata", paste0("ramsar_site_data_100m_", c_name))
  if (dir.exists(data_dir)) {
    countries <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
    for (country in countries) {
      country_to_cont[[country]] <- c_name
    }
  }
}

get_site_continent <- function(country_name) {
  cont <- country_to_cont[[country_name]]
  if (!is.null(cont)) return(cont)
  # Fallback: case-insensitive match
  names_lower <- tolower(names(country_to_cont))
  idx <- match(tolower(country_name), names_lower)
  if (!is.na(idx)) return(country_to_cont[[idx]])
  return("unknown")
}

wkt_files <- list.files(ramsar_wkt_base_dir, pattern = "\\.wkt$", recursive = TRUE, full.names = TRUE)
target_wkt_files <- c()
for (fid in target_ramsarids) {
  match_pat <- paste0("^site_", fid, "_")
  matched <- wkt_files[grepl(match_pat, basename(wkt_files))]
  target_wkt_files <- c(target_wkt_files, matched)
}
target_wkt_files <- unique(target_wkt_files)
cat("Found", length(target_wkt_files), "matching WKT files for target sites.\n")

if (length(target_wkt_files) == 0) {
  stop("No target WKT files found. Ensure step1 has run successfully.")
}

ramsar_sf_list <- lapply(target_wkt_files, function(file_path) {
  tryCatch({
    wkt_string <- readLines(file_path, warn = FALSE)
    geom <- st_as_sfc(wkt_string, crs = 4326)
    filename <- basename(file_path)
    ramsar_site_id <- tools::file_path_sans_ext(filename)
    country_name <- basename(dirname(file_path))
    site_continent <- get_site_continent(country_name)
    st_as_sf(data.frame(
      ramsar_site_id = ramsar_site_id,
      site_name = ramsar_site_id,
      country_en = country_name,
      continent = site_continent,
      source_wkt_file = filename
    ), geom = geom)
  }, error = function(e) {
    cat("  [Warning] Failed to process WKT file:", file_path, "-", e$message, "\n")
    return(NULL)
  })
})

ramsar_sites_sf <- do.call(rbind, Filter(Negate(is.null), ramsar_sf_list))
ramsar_sites_sf <- st_make_valid(ramsar_sites_sf)
cat("Target site geometries loaded and validated.\n")

# --- Step 3.5: Extract target 100km MGRS square IDs by continent ---
cat("\n=== Extracting target 100km MGRS square IDs by continent ===\n")
target_squares_by_cont <- list()
for (c_name in continents_list) {
  target_squares_by_cont[[c_name]] <- c()
}

original_s2_setting <- sf_use_s2()
sf_use_s2(FALSE)

for (i in 1:nrow(ramsar_sites_sf)) {
  geom <- st_geometry(ramsar_sites_sf[i, ])
  bbox <- st_bbox(geom)
  site_cont <- ramsar_sites_sf$continent[i]
  
  # Generate grid of points with spacing of ~5.5km (0.05 degrees)
  x_step <- max(0.01, (bbox["xmax"] - bbox["xmin"]) / 20)
  y_step <- max(0.01, (bbox["ymax"] - bbox["ymin"]) / 20)
  
  xs <- seq(bbox["xmin"], bbox["xmax"], by = x_step)
  ys <- seq(bbox["ymin"], bbox["ymax"], by = y_step)
  grid_pts <- expand.grid(lng = xs, lat = ys)
  
  # Vertices
  vertices <- st_coordinates(geom)
  pts_all <- rbind(
    grid_pts,
    data.frame(lng = vertices[, "X"], lat = vertices[, "Y"])
  )
  
  pts_all$lat_round <- round(pts_all$lat, 4)
  pts_all$lng_round <- round(pts_all$lng, 4)
  unique_pts <- unique(pts_all[, c("lat_round", "lng_round")])
  
  sqs <- tryCatch({
    mapply(mgrs::latlng_to_mgrs, unique_pts$lat_round, unique_pts$lng_round, MoreArgs = list(precision = 0))
  }, error = function(e) NULL)
  
  if (!is.null(sqs)) {
    sqs <- unique(na.omit(sqs))
    if (site_cont %in% continents_list) {
      target_squares_by_cont[[site_cont]] <- unique(c(target_squares_by_cont[[site_cont]], sqs))
    } else {
      for (c_name in continents_list) {
        target_squares_by_cont[[c_name]] <- unique(c(target_squares_by_cont[[c_name]], sqs))
      }
    }
  }
}
sf_use_s2(original_s2_setting)

for (c_name in continents_list) {
  cat(sprintf("  %s: %d unique MGRS squares\n", c_name, length(target_squares_by_cont[[c_name]])))
}

# --- Step 4: Process continental cubes ---
for (cont in continents) {
  cube_path <- file.path("inst/extdata/continental_gbif_data", paste0(cont, "_100m_mgrs_all.csv"))
  if (!file.exists(cube_path)) {
    cat("\nContinental cube not found for:", cont, "- skipping.\n")
    next
  }
  
  output_dir <- file.path("inst/extdata", paste0("ramsar_site_data_100m_", cont))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat(sprintf("\n=== Processing Continent: %s ===\n", toupper(cont)))
  cat("Reading header of:", cube_path, "\n")
  header <- names(fread(cube_path, nrows = 0))
  if (!(mgrs_column_name %in% header)) {
    cat("Error: MGRS column not found in", cube_path, "\n")
    next
  }
  
  total_rows_processed <- resume_rows_by_cont[[cont]]
  chunk_num <- floor(total_rows_processed / chunk_size) + 1
  
  repeat {
    withr::local_options(list(scipen = 999))
    cat(sprintf("\n- Continent %s, Chunk %d (start line %s)...\n", cont, chunk_num, format(total_rows_processed + 2, scientific=FALSE)))
    
    data_chunk <- tryCatch({
      if (total_rows_processed == 0) {
        fread(cube_path, nrows = chunk_size, encoding = "UTF-8", showProgress = FALSE)
      } else {
        # Pure R fread with skip to avoid external shell tail/head utilities on Windows
        fread(cube_path, skip = total_rows_processed + 1, nrows = chunk_size, header = FALSE, col.names = header, encoding = "UTF-8", showProgress = FALSE)
      }
    }, error = function(e) {
      cat("  [Error] Failed to read chunk:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(data_chunk) || nrow(data_chunk) == 0) {
      cat("- Reached end of file.\n")
      break
    }
    
    unique_mgrs_codes <- unique(data_chunk[[mgrs_column_name]])
    # Sanitize unique MGRS codes to prevent multibyte string parsing errors
    unique_mgrs_codes <- iconv(unique_mgrs_codes, to = "UTF-8", sub = "")
    unique_prefixes <- substr(unique_mgrs_codes, 1, 5)
    matching_idx <- unique_prefixes %in% target_squares_by_cont[[cont]]
    filtered_mgrs_codes <- unique_mgrs_codes[matching_idx]
    
    cat(sprintf("  Pre-filtering chunk: keeping %d of %d unique MGRS cells (%.2f%%) matching target squares.\n", 
                length(filtered_mgrs_codes), length(unique_mgrs_codes),
                (length(filtered_mgrs_codes)/length(unique_mgrs_codes))*100))
    
    if (length(filtered_mgrs_codes) == 0) {
      cat("  No matching cells in this chunk. Skipping.\n")
      total_rows_processed <- total_rows_processed + nrow(data_chunk)
      chunk_num <- chunk_num + 1
      next
    }
    
    cat(sprintf("  Converting %d unique MGRS cells to lat/lon...\n", length(filtered_mgrs_codes)))
    latlon_coords <- tryCatch({
      mgrs_to_latlng(filtered_mgrs_codes)
    }, error = function(e) {
      cat("  [Warning] MGRS conversion error:", e$message, "\n")
      return(data.frame(mgrs=character(), lat=numeric(), lon=numeric()))
    })
    
    valid_coords <- na.omit(latlon_coords)
    if (nrow(valid_coords) == 0) {
      cat("  No valid coordinates in chunk. Skipping.\n")
      total_rows_processed <- total_rows_processed + nrow(data_chunk)
      chunk_num <- chunk_num + 1
      next
    }
    
    points_sf <- st_as_sf(valid_coords, coords = c("lng", "lat"), crs = 4326)
    
    cat("  Intersecting with target Ramsar site geometries...\n")
    original_s2_setting <- sf_use_s2()
    sf_use_s2(FALSE)
    points_in_ramsar <- st_join(points_sf, ramsar_sites_sf, join = st_intersects, left = FALSE)
    sf_use_s2(original_s2_setting)
    
    if (nrow(points_in_ramsar) == 0) {
      cat("  No records matched target sites in this chunk.\n")
    } else {
      cat(sprintf("  Matched %d records to target sites. Writing files...\n", nrow(points_in_ramsar)))
      setnames(points_in_ramsar, "mgrs", mgrs_column_name)
      
      results_to_write_dt <- merge(
        as.data.table(points_in_ramsar),
        as.data.table(data_chunk),
        by = mgrs_column_name
      )
      
      split_by_site_country <- split(results_to_write_dt, by = c("ramsar_site_id", "country_en"))
      
      for (site_country_key in names(split_by_site_country)) {
        site_data <- split_by_site_country[[site_country_key]]
        key_parts <- strsplit(site_country_key, "\\.")[[1]]
        site_id <- key_parts[1]
        country_name <- key_parts[2]
        
        country_output_dir <- file.path(output_dir, country_name)
        if (!dir.exists(country_output_dir)) {
          dir.create(country_output_dir, recursive = TRUE)
        }
        
        output_filename <- file.path(country_output_dir, paste0(site_id, "_data.csv"))
        file_exists_on_disk <- file.exists(output_filename)
        
        fwrite(site_data,
               file = output_filename,
               append = file_exists_on_disk,
               col.names = !file_exists_on_disk)
        
        cat(sprintf("    %s -> wrote %d records (append = %s)\n", basename(output_filename), nrow(site_data), file_exists_on_disk))
      }
    }
    
    total_rows_processed <- total_rows_processed + nrow(data_chunk)
    chunk_num <- chunk_num + 1
  }
}

cat("\n=== Targeted Cube Splitting Process Complete for Remaining Continents ===\n")
