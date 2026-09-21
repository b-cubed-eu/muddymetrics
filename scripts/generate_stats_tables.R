#!/usr/bin/env Rscript
library(dplyr)

p <- read.csv("output/data_sufficiency/subgroup_analysis/sites_passing_all_precision.csv")
full <- read.csv("output/data_sufficiency/subgroup_analysis/all_tests_results_with_precision.csv")

# --- TABLE 1: Passing sites by continent ---
cat("=== TABLE 1: PASSING SITES BY CONTINENT ===\n")
by_cont <- p %>% group_by(continent) %>%
  summarise(combos = n(), unique_sites = n_distinct(site_id), .groups = "drop") %>%
  arrange(desc(unique_sites))
print(by_cont)

# --- TABLE 2: Total sites evaluated per continent ---
cat("\n=== TABLE 2: TOTAL SITES EVALUATED PER CONTINENT ===\n")
total_by_cont <- full %>% group_by(continent) %>%
  summarise(total_sites = n_distinct(site_id), total_combos = n(), .groups = "drop")
print(total_by_cont)

# --- TABLE 3: Invisibility metric ---
cat("\n=== TABLE 3: INVISIBILITY METRIC BY CONTINENT ===\n")
inv <- total_by_cont %>%
  left_join(by_cont %>% select(continent, unique_sites), by = "continent") %>%
  mutate(
    unique_sites = ifelse(is.na(unique_sites), 0, unique_sites),
    invisible_pct = round((1 - unique_sites / total_sites) * 100, 1)
  )
print(inv)
cat("\nGLOBAL:", sum(inv$total_sites), "sites evaluated,", sum(inv$unique_sites), "pass all criteria\n")
cat("Global Invisibility:", round((1 - sum(inv$unique_sites) / sum(inv$total_sites)) * 100, 1), "%\n")

# --- TABLE 4: By taxonomic group ---
cat("\n=== TABLE 4: BY TAXONOMIC GROUP ===\n")
print(p %>% group_by(subset) %>% summarise(combos = n(), .groups = "drop") %>% arrange(desc(combos)))

# --- TABLE 5: Precision ratio stats ---
cat("\n=== TABLE 5: PRECISION RATIO STATS (ALL EVALUATED) ===\n")
cat("Mean:", round(mean(full$precision_ratio, na.rm = TRUE), 4), "\n")
cat("Median:", round(median(full$precision_ratio, na.rm = TRUE), 4), "\n")
cat("% with NA precision:", round(mean(is.na(full$precision_ratio)) * 100, 1), "\n")

# --- TABLE 6: Individual test pass rates (global) ---
cat("\n=== TABLE 6: INDIVIDUAL TEST PASS RATES (GLOBAL) ===\n")
cat("Density:   ", round(mean(full$pass_density, na.rm = TRUE) * 100, 1), "%\n")
cat("Chao2:     ", round(mean(full$pass_chao2, na.rm = TRUE) * 100, 1), "%\n")
cat("SAC Slope: ", round(mean(full$pass_sac_slope, na.rm = TRUE) * 100, 1), "%\n")
cat("Temporal:  ", round(mean(full$pass_temporal, na.rm = TRUE) * 100, 1), "%\n")
cat("Precision: ", round(mean(full$pass_precision, na.rm = TRUE) * 100, 1), "%\n")
cat("ALL:       ", round(mean(full$pass_all, na.rm = TRUE) * 100, 1), "%\n")

# --- TABLE 7: Top 15 passing sites with details ---
cat("\n=== TABLE 7: TOP PASSING SITES (BY # SUBGROUPS PASSING) ===\n")
top_sites <- p %>%
  group_by(continent, country, site_id) %>%
  summarise(
    subgroups_passing = n(),
    subsets = paste(subset, collapse = ", "),
    avg_density = round(mean(density, na.rm = TRUE), 2),
    avg_chao2 = round(mean(chao2, na.rm = TRUE), 2),
    avg_precision = round(mean(precision_ratio, na.rm = TRUE), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(subgroups_passing))
print(top_sites, n = 20)

# --- TABLE 8: Bottleneck analysis ---
cat("\n=== TABLE 8: BOTTLENECK ANALYSIS ===\n")
cat("Which test eliminates the most site-subgroups?\n")
# Start from all 6675, progressively filter
n_total <- nrow(full)
n_pass_density <- sum(full$pass_density, na.rm = TRUE)
n_pass_chao2 <- sum(full$pass_density & full$pass_chao2, na.rm = TRUE)
n_pass_sac <- sum(full$pass_density & full$pass_chao2 & full$pass_sac_slope, na.rm = TRUE)
n_pass_temp <- sum(full$pass_density & full$pass_chao2 & full$pass_sac_slope & full$pass_temporal, na.rm = TRUE)
n_pass_prec <- sum(full$pass_all, na.rm = TRUE)

cat(sprintf("Start:                %d combos\n", n_total))
cat(sprintf("After Density:        %d (-%d, %.1f%% eliminated)\n", n_pass_density, n_total - n_pass_density, (n_total - n_pass_density)/n_total*100))
cat(sprintf("After + Chao2:        %d (-%d, %.1f%% eliminated)\n", n_pass_chao2, n_pass_density - n_pass_chao2, (n_pass_density - n_pass_chao2)/n_pass_density*100))
cat(sprintf("After + SAC Slope:    %d (-%d, %.1f%% eliminated)\n", n_pass_sac, n_pass_chao2 - n_pass_sac, (n_pass_chao2 - n_pass_sac)/n_pass_chao2*100))
cat(sprintf("After + Temporal:     %d (-%d, %.1f%% eliminated)\n", n_pass_temp, n_pass_sac - n_pass_temp, (n_pass_sac - n_pass_temp)/n_pass_sac*100))
cat(sprintf("After + Precision:    %d (-%d, %.1f%% eliminated)\n", n_pass_prec, n_pass_temp - n_pass_prec, (n_pass_temp - n_pass_prec)/n_pass_temp*100))

# --- TABLE 9: Countries with most passing sites ---
cat("\n=== TABLE 9: TOP COUNTRIES BY PASSING SITES ===\n")
country_stats <- p %>%
  group_by(country) %>%
  summarise(n_sites = n_distinct(site_id), n_combos = n(), .groups = "drop") %>%
  arrange(desc(n_sites))
print(country_stats, n = 20)
