library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(impIndicator)
library(gridExtra)

# Sizing & Theme Constants
ACCENT_BLUE <- "#1A73E8"
ACCENT_RED <- "#D93025"
ACCENT_GREEN <- "#188038"
NEUTRAL_GREY <- "#70757A"
TITLE_COLOR <- "#202124"
TEXT_COLOR <- "#3C4043"

base_theme <- theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 12, color = TITLE_COLOR, hjust = 0.5),
    plot.subtitle = element_text(size = 9, color = TEXT_COLOR, hjust = 0.5),
    axis.title = element_text(face = "bold", size = 9, color = TEXT_COLOR),
    axis.text = element_text(size = 8, color = TEXT_COLOR),
    legend.title = element_text(face = "bold", size = 8, color = TEXT_COLOR),
    legend.text = element_text(size = 8, color = TEXT_COLOR),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

gidias_file <- "inst/extdata/GIDIAS/GIDIAS_machine_read.csv"
out_dir <- "report_figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Loading GIDIAS database...\n")
gidias <- read_csv(gidias_file, show_col_types = FALSE)

# Define the sites to analyze
sites <- list(
  france = list(
    id = 786,
    name = "La Petite Camargue (France)",
    file = "temp_worker_full_site_786_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/france/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/france/speciesprofile.txt"
  ),
  uk = list(
    id = 67,
    name = "Site 67 (United Kingdom)",
    file = "temp_worker_full_site_67_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/united_kingdom/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/united_kingdom/speciesprofile.txt"
  ),
  botswana = list(
    id = 879,
    name = "Site 879 (Botswana)",
    file = "temp_worker_full_site_879_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/Botswana/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/Botswana/speciesprofile.txt"
  )
)

combined_results <- list()

for (skey in names(sites)) {
  s <- sites[[skey]]
  cat(sprintf("\n--- Processing site %d: %s ---\n", s$id, s$name))
  
  # 1. Load GRIIS Data
  cat("Loading GRIIS data...\n")
  taxon_df <- read.delim(s$griis_taxon, stringsAsFactors = FALSE, sep = "\t", quote = "")
  profile_df <- read.delim(s$griis_profile, stringsAsFactors = FALSE, sep = "\t", quote = "")
  
  taxon_df <- taxon_df %>%
    mutate(clean_name = str_extract(scientificName, "^[^ ]+ [^ ]+"))
  
  alien_registry <- taxon_df %>%
    left_join(profile_df, by = "id") %>%
    select(clean_name, isInvasive) %>%
    distinct(clean_name, .keep_all = TRUE)
  
  # 2. Load Site Occurrence Cube
  cat("Loading site occurrence cube...\n")
  site_cube <- read_csv(s$file, show_col_types = FALSE)
  
  # Intersect with GRIIS first
  site_mapped <- site_cube %>%
    left_join(alien_registry, by = c("species" = "clean_name")) %>%
    mutate(
      is_alien = !is.na(isInvasive),
      is_invasive = !is.na(isInvasive) & isInvasive == "Invasive"
    )
  
  alien_species_list <- unique(site_mapped$species[site_mapped$is_alien])
  invasive_species_list <- unique(site_mapped$species[site_mapped$is_invasive])
  
  cat(sprintf("Site species: %d | Alien: %d | Invasive: %d\n", 
              length(unique(site_cube$species)), length(alien_species_list), length(invasive_species_list)))
  
  # 3. Match to GIDIAS EICAT Data
  eicat_data <- gidias %>%
    filter(Verified.Name.GBIF.Taxon %in% alien_species_list) %>%
    filter(!is.na(magnitude.Nature)) %>%
    mutate(
      scientific_name = Verified.Name.GBIF.Taxon,
      impact_category = case_when(
        magnitude.Nature == 0 ~ "MC",
        magnitude.Nature == 1 ~ "MN",
        magnitude.Nature == 2 ~ "MO",
        magnitude.Nature == 3 ~ "MR",
        magnitude.Nature == 4 ~ "MV",
        TRUE ~ "MC"
      ),
      impact_mechanism = coalesce(mechanism.Nature.clean, "Unknown")
    ) %>%
    select(scientific_name, impact_category, impact_mechanism) %>%
    distinct()
  
  cat(sprintf("Matched to GIDIAS: %d alien species\n", length(unique(eicat_data$scientific_name))))
  
  # Run impIndicator overall trend
  site_cube_imp <- site_cube %>%
    rename(
      scientificName = species,
      cellCode = mgrscellcode,
      taxonKey = specieskey
    )
  
  cat("Computing standardized EICAT impact indicator...\n")
  res_imp <- compute_impact_indicator(
    cube = site_cube_imp,
    impact_data = eicat_data,
    method = "mean_cum",
    ci_type = "none"
  )
  impact_df <- res_imp$impact %>%
    rename(impact_indicator = diversity_val)
  
  # Match and aggregate annual stats
  annual_stats <- site_mapped %>%
    group_by(year) %>%
    summarise(
      total_richness = n_distinct(species),
      alien_richness = n_distinct(species[is_alien]),
      invasive_richness = n_distinct(species[is_invasive]),
      total_occurrences = sum(occurrences),
      .groups = "drop"
    ) %>%
    left_join(impact_df, by = "year") %>%
    mutate(
      site_id = s$id,
      site_name = s$name
    )
  
  # Calculate correlation between EICAT impact indicator and total occurrences
  clean_stats <- annual_stats %>% filter(!is.na(impact_indicator), !is.na(total_occurrences))
  if (nrow(clean_stats) > 2) {
    corr_val <- cor(clean_stats$impact_indicator, clean_stats$total_occurrences)
    cat(sprintf("Correlation (EICAT Impact vs. Total Occurrences): %.4f\n", corr_val))
  } else {
    cat("Insufficient years for correlation calculation.\n")
  }
  
  write_csv(annual_stats, file.path(out_dir, sprintf("site_%d_invasive_impact.csv", s$id)))
  
  combined_results[[skey]] <- annual_stats
}

# Combine all results for comparison plotting
all_sites_df <- bind_rows(combined_results)

cat("\nGenerating comparative plot...\n")
# We'll create a multi-panel plot to show EICAT impact vs occurrences for all three sites
p_comp <- ggplot(all_sites_df, aes(x = year)) +
  geom_line(aes(y = impact_indicator, color = "EICAT Impact Indicator"), linewidth = 1.2) +
  geom_point(aes(y = impact_indicator, color = "EICAT Impact Indicator"), size = 2) +
  # Add sampling effort (occurrences) on a secondary scale
  # We scale the occurrences for each site to fit the 0-max(impact_indicator) scale
  geom_line(aes(y = total_occurrences / 100000 * 0.3, color = "Total Occurrences (scaled)"), linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~site_name, ncol = 1, scales = "free_y") +
  scale_color_manual(
    name = NULL,
    values = c(
      "EICAT Impact Indicator" = ACCENT_RED,
      "Total Occurrences (scaled)" = NEUTRAL_GREY
    )
  ) +
  labs(
    title = "Comparison of EICAT Invasive Impact & Sampling Effort across Sites",
    subtitle = "Comparing standardized cumulative EICAT impact indicator trends to raw occurrence records",
    x = "Year",
    y = "Impact Indicator Value"
  ) +
  base_theme +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "comparative_sites_invasive_impact.png"), plot = p_comp, width = 10.0, height = 8.5, dpi = 300)
ggsave(file.path(out_dir, "comparative_sites_invasive_impact.svg"), plot = p_comp, width = 10.0, height = 8.5)

cat("\nAll sites processed successfully! Comparison figures saved.\n")
