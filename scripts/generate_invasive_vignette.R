library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(impIndicator)

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

# Paths
griis_taxon_file <- "inst/extdata/griis_checklists/france/taxon.txt"
griis_profile_file <- "inst/extdata/griis_checklists/france/speciesprofile.txt"
site_data_file <- "temp_worker_full_site_786_full.csv"
gidias_file <- "inst/extdata/GIDIAS/GIDIAS_machine_read.csv"
out_dir <- "report_figures"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 1. Load GRIIS Data
cat("Loading GRIIS data...\n")
taxon_df <- read.delim(griis_taxon_file, stringsAsFactors = FALSE, sep = "\t", quote = "")
profile_df <- read.delim(griis_profile_file, stringsAsFactors = FALSE, sep = "\t", quote = "")

taxon_df <- taxon_df %>%
  mutate(clean_name = str_extract(scientificName, "^[^ ]+ [^ ]+"))

alien_registry <- taxon_df %>%
  left_join(profile_df, by = "id") %>%
  select(clean_name, isInvasive) %>%
  distinct(clean_name, .keep_all = TRUE)

# 2. Load Site Occurrence Cube
cat("Loading site 786 occurrence cube...\n")
site_cube <- read_csv(site_data_file, show_col_types = FALSE)

# 3. Load GIDIAS EICAT Data
cat("Loading GIDIAS data...\n")
gidias <- read_csv(gidias_file, show_col_types = FALSE)

# Build list of alien species at site 786 based on GRIIS France
site_mapped <- site_cube %>%
  left_join(alien_registry, by = c("species" = "clean_name")) %>%
  mutate(
    is_alien = !is.na(isInvasive),
    is_invasive = !is.na(isInvasive) & isInvasive == "Invasive"
  )

alien_species_list <- unique(site_mapped$species[site_mapped$is_alien])

# Filter GIDIAS EICAT data only for site_786 alien species
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

# 4. Generate Figure 7: Raw Richness vs. impIndicator EICAT trend
cat("Generating Figure 7...\n")
# Match and aggregate annual stats
site_mapped <- site_cube %>%
  left_join(alien_registry, by = c("species" = "clean_name")) %>%
  mutate(
    is_alien = !is.na(isInvasive),
    is_invasive = !is.na(isInvasive) & isInvasive == "Invasive"
  )

annual_stats <- site_mapped %>%
  group_by(year) %>%
  summarise(
    total_richness = n_distinct(species),
    alien_richness = n_distinct(species[is_alien]),
    invasive_richness = n_distinct(species[is_invasive]),
    total_occurrences = sum(occurrences),
    .groups = "drop"
  ) %>%
  left_join(impact_df, by = "year")

write_csv(annual_stats, file.path(out_dir, "france_site_786_invasive_impact.csv"))

# Panel A: Richness trends over time
p7a <- ggplot(annual_stats, aes(x = year)) +
  geom_line(aes(y = total_richness / 50, color = "Total Richness (scaled 1/50)"), linewidth = 1.0, linetype = "dashed") +
  geom_line(aes(y = alien_richness, color = "Alien Richness"), linewidth = 1.2) +
  geom_point(aes(y = alien_richness, color = "Alien Richness"), size = 2) +
  geom_line(aes(y = invasive_richness, color = "Invasive Richness"), linewidth = 1.2) +
  geom_point(aes(y = invasive_richness, color = "Invasive Richness"), size = 2) +
  scale_y_continuous(
    name = "Alien / Invasive Richness",
    sec.axis = sec_axis(~.*50, name = "Total Species Richness")
  ) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Alien Richness" = ACCENT_BLUE,
      "Invasive Richness" = ACCENT_RED,
      "Total Richness (scaled 1/50)" = NEUTRAL_GREY
    )
  ) +
  labs(
    title = "Panel A: Recorded Richness & Sampling Effort Trends",
    subtitle = "Raw species counts rise with sampling volume",
    x = "Year",
    y = "Species Richness"
  ) +
  base_theme +
  theme(legend.position = "top")

# Panel B: EICAT standardized impact indicator
p7b <- ggplot(annual_stats, aes(x = year, y = impact_indicator)) +
  geom_line(color = ACCENT_RED, linewidth = 1.2) +
  geom_point(color = ACCENT_RED, size = 2) +
  labs(
    title = "Panel B: Standardised EICAT Invasive Impact Trend",
    subtitle = "impIndicator cumulative impact value (method 'mean_cum')",
    x = "Year",
    y = "Impact Indicator Value"
  ) +
  base_theme

# Combine into a single vertical multi-panel plot
library(gridExtra)
p7 <- grid.arrange(p7a, p7b, ncol = 1, heights = c(1, 1))

ggsave(file.path(out_dir, "france_site_786_invasive_trend.png"), plot = p7, width = 9.0, height = 7.5, dpi = 300)
ggsave(file.path(out_dir, "france_site_786_invasive_trend.svg"), plot = p7, width = 9.0, height = 7.5)

# 5. Generate Figure 8: Horizontal bar plot of top 10 invasive species color-coded by EICAT
cat("Generating Figure 8...\n")
top_10 <- site_mapped %>%
  filter(is_invasive) %>%
  group_by(species) %>%
  summarise(total_occ = sum(occurrences), .groups = "drop") %>%
  arrange(desc(total_occ)) %>%
  head(10)

# Fetch maximum EICAT category for each of the top 10 species from eicat_data
top_10_eicat <- top_10 %>%
  left_join(
    eicat_data %>%
      group_by(scientific_name) %>%
      summarise(
        impact_category = case_when(
          "MV" %in% impact_category ~ "Massive (MV)",
          "MR" %in% impact_category ~ "Major (MR)",
          "MO" %in% impact_category ~ "Moderate (MO)",
          "MN" %in% impact_category ~ "Minor (MN)",
          "MC" %in% impact_category ~ "Minimal (MC)",
          TRUE ~ "Not assessed"
        ),
        .groups = "drop"
      ),
    by = c("species" = "scientific_name")
  ) %>%
  mutate(
    impact_category = factor(coalesce(impact_category, "Not assessed"), 
                             levels = c("Not assessed", "Minimal (MC)", "Minor (MN)", "Moderate (MO)", "Major (MR)", "Massive (MV)"))
  )

# Add display names matching Table 8 in the report
display_names <- c(
  "Cortaderia selloana" = "Cortaderia selloana\n(pampas grass)",
  "Amorpha fruticosa" = "Amorpha fruticosa\n(false indigo-bush)",
  "Myocastor coypus" = "Myocastor coypus\n(coypu)",
  "Pseudorasbora parva" = "Pseudorasbora parva\n(topmouth gudgeon)",
  "Baccharis halimifolia" = "Baccharis halimifolia\n(saltbush)",
  "Ludwigia peploides" = "Ludwigia peploides\n(water primrose)",
  "Procambarus clarkii" = "Procambarus clarkii\n(red swamp crayfish)",
  "Robinia pseudoacacia" = "Robinia pseudoacacia\n(black locust)",
  "Trachemys scripta" = "Trachemys scripta\n(pond slider)",
  "Cyprinus carpio" = "Cyprinus carpio\n(common carp)"
)

top_10_eicat <- top_10_eicat %>%
  mutate(
    display_name = display_names[species],
    display_name = reorder(display_name, total_occ)
  )

p8 <- ggplot(top_10_eicat, aes(x = total_occ, y = display_name, fill = impact_category)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.2) +
  geom_text(aes(label = scales::comma(total_occ)), hjust = -0.15, size = 3, fontface = "bold", color = TEXT_COLOR) +
  scale_fill_manual(
    name = "Maximum EICAT Impact Category",
    values = c(
      "Not assessed" = "#E0E0E0",
      "Minimal (MC)" = "#9AA0A6",
      "Minor (MN)" = ACCENT_BLUE,
      "Moderate (MO)" = "#FBBC04",
      "Major (MR)" = "#E37400",
      "Massive (MV)" = ACCENT_RED
    ),
    drop = FALSE
  ) +
  scale_x_continuous(
    name = "Total Occurrences (GBIF, 1990-2025)",
    limits = c(0, max(top_10_eicat$total_occ) * 1.15),
    labels = scales::comma
  ) +
  labs(
    title = "Most Frequently Recorded Invasive Taxa at La Petite Camargue",
    subtitle = "Top 10 invasive species by occurrences, color-coded by maximum GIDIAS EICAT impact category",
    y = ""
  ) +
  base_theme +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    plot.margin = margin(15, 30, 15, 15)
  )

ggsave(file.path(out_dir, "France_site_786_full_top_aliens.png"), plot = p8, width = 10.0, height = 6.0, dpi = 300)
ggsave(file.path(out_dir, "France_site_786_full_top_aliens.svg"), plot = p8, width = 10.0, height = 6.0)

cat("Done! Figure 7 and Figure 8 generated in report_figures/.\n")
