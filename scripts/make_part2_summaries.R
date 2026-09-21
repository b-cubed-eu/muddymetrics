#!/usr/bin/env Rscript
library(dplyr)
library(purrr)
library(tools)

# Load passing combos and MAGV summary
passing_file <- "output/data_sufficiency/subgroup_analysis/sites_passing_corrected.csv"
magv_file <- "output/magv_summary.csv"

if (!file.exists(passing_file) || !file.exists(magv_file)) {
  stop("Input files sites_passing_corrected.csv or magv_summary.csv not found!")
}

passing_combos <- read.csv(passing_file, stringsAsFactors = FALSE)
magv_summary <- read.csv(magv_file, stringsAsFactors = FALSE)

# Helper to find data frame in nested lists
find_df <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x)) {
    for (n in names(x)) {
      res <- find_df(x[[n]])
      if (!is.null(res)) return(res)
    }
  }
  return(NULL)
}

indicator_records <- list()

cat("Processing full analysis RDS files...\n")
for (i in seq_len(nrow(passing_combos))) {
  row <- passing_combos[i, ]
  site_id <- row$site_id
  sub <- row$subset
  cont <- row$continent
  country <- row$country
  
  # Find stable window info
  magv_row <- magv_summary %>%
    filter(site_id == !!site_id & subset == !!sub)
  
  has_stable <- FALSE
  stable_start <- NA_integer_
  stable_end <- NA_integer_
  if (nrow(magv_row) > 0) {
    has_stable <- magv_row$has_usable_window[1]
    stable_start <- magv_row$stable_run_start[1]
    stable_end <- magv_row$stable_run_end[1]
  }
  
  site_dir <- file.path("output/full_analysis", site_id)
  if (!dir.exists(site_dir)) next
  
  rds_files <- list.files(site_dir, pattern = paste0("^", sub, "_.*_ts\\.rds$"), full.names = TRUE)
  
  for (rf in rds_files) {
    # Extract indicator name from filename
    # E.g., aves_obs_richness_ts.rds -> obs_richness
    base_name <- file_path_sans_ext(basename(rf))
    parts <- strsplit(base_name, "_")[[1]]
    indicator_id <- paste(parts[2:(length(parts)-1)], collapse = "_")
    
    data_obj <- tryCatch(readRDS(rf), error = function(e) NULL)
    if (is.null(data_obj)) next
    
    df <- find_df(data_obj)
    if (is.null(df) || nrow(df) == 0) next
    
    # Identify value column
    val_cols <- c("diversity_val", "value", "richness", "diversity", "evenness", "completeness", 
                  "density", "hill0", "hill1", "hill2", "pielou", "williams")
    existing_cols <- intersect(names(df), val_cols)
    if (length(existing_cols) == 0) next
    val_col <- existing_cols[1]
    
    # Sort by year
    df <- df %>% arrange(year)
    
    n_yrs <- nrow(df)
    first_yr <- min(df$year)
    last_yr <- max(df$year)
    mean_val <- mean(df[[val_col]], na.rm = TRUE)
    start_val <- df[[val_col]][1]
    end_val <- df[[val_col]][n_yrs]
    
    # Linear model over full time series
    trend_slope <- NA_real_
    trend_p <- NA_real_
    if (n_yrs >= 3) {
      lm_fit <- tryCatch(lm(df[[val_col]] ~ df$year), error = function(e) NULL)
      if (!is.null(lm_fit)) {
        summary_lm <- summary(lm_fit)
        if (nrow(summary_lm$coefficients) >= 2) {
          trend_slope <- summary_lm$coefficients[2, 1]
          trend_p <- summary_lm$coefficients[2, 4]
        }
      }
    }
    
    # CIs
    ci_lower_mean <- NA_real_
    ci_upper_mean <- NA_real_
    if ("ll" %in% names(df)) {
      ci_lower_mean <- mean(df$ll, na.rm = TRUE)
    }
    if ("ul" %in% names(df)) {
      ci_upper_mean <- mean(df$ul, na.rm = TRUE)
    }
    
    # Stable window trend
    trend_slope_stable <- NA_real_
    first_yr_stable <- NA_integer_
    last_yr_stable <- NA_integer_
    
    if (has_stable && !is.na(stable_start) && !is.na(stable_end)) {
      df_stable <- df %>% filter(year >= stable_start & year <= stable_end)
      n_stable_yrs <- nrow(df_stable)
      if (n_stable_yrs >= 1) {
        first_yr_stable <- min(df_stable$year)
        last_yr_stable <- max(df_stable$year)
      }
      if (n_stable_yrs >= 3) {
        lm_fit_stable <- tryCatch(lm(df_stable[[val_col]] ~ df_stable$year), error = function(e) NULL)
        if (!is.null(lm_fit_stable)) {
          summary_lm_stable <- summary(lm_fit_stable)
          if (nrow(summary_lm_stable$coefficients) >= 2) {
            trend_slope_stable <- summary_lm_stable$coefficients[2, 1]
          }
        }
      }
    }
    
    indicator_records[[length(indicator_records) + 1]] <- data.frame(
      site_id = site_id,
      subset = sub,
      continent = cont,
      country = country,
      indicator = indicator_id,
      n_years = n_yrs,
      first_year = first_yr,
      last_year = last_yr,
      mean_value = mean_val,
      start_value = start_val,
      end_value = end_val,
      trend_slope = trend_slope,
      trend_p_value = trend_p,
      ci_lower_mean = ci_lower_mean,
      ci_upper_mean = ci_upper_mean,
      trend_slope_stable = trend_slope_stable,
      first_year_stable = first_yr_stable,
      last_year_stable = last_yr_stable,
      stringsAsFactors = FALSE
    )
  }
}

if (length(indicator_records) > 0) {
  indicator_summary_df <- do.call(rbind, indicator_records)
  write.csv(indicator_summary_df, "output/part2_indicator_summary.csv", row.names = FALSE)
  cat(sprintf("✓ Saved output/part2_indicator_summary.csv with %d records.\n", nrow(indicator_summary_df)))
  
  # Generate trend direction summary
  # Classify: Positive (slope > 0, p <= 0.1), Negative (slope < 0, p <= 0.1), Neutral-Stable (p > 0.1 or slope == 0)
  trend_direction_df <- indicator_summary_df %>%
    mutate(
      direction = case_when(
        is.na(trend_slope) ~ "Neutral",
        is.na(trend_p_value) ~ "Neutral",
        trend_p_value <= 0.1 & trend_slope > 0 ~ "Positive",
        trend_p_value <= 0.1 & trend_slope < 0 ~ "Negative",
        TRUE ~ "Neutral"
      )
    ) %>%
    group_by(continent, subset, indicator) %>%
    summarise(
      n_combos = n(),
      n_positive = sum(direction == "Positive"),
      n_negative = sum(direction == "Negative"),
      n_neutral = sum(direction == "Neutral"),
      .groups = "drop"
    )
  
  write.csv(trend_direction_df, "output/part2_trend_direction.csv", row.names = FALSE)
  cat(sprintf("✓ Saved output/part2_trend_direction.csv with %d rows.\n", nrow(trend_direction_df)))
  
} else {
  cat("No indicator records found to summarize.\n")
}
