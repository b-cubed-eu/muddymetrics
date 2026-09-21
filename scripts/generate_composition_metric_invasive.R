# ---------------------------------------------------------------------------
# Ratio-based composition metric for invasive impact (Part 3 extension)
#
# Motivation: both rarefactions (size-based and coverage-based) showed that
# impIndicator's `mean_cum` output is close to a monotone function of sample
# size - it sums per-cell impact means across occupied cells, so it grows with
# the number of cells a sample touches. See coverage_standardisation_findings.md.
#
# A ratio has no built-in sample-size dependence: it asks what FRACTION of the
# alien records are high-impact, not how many there are. No rarefaction needed.
#
#   prop_records_major_plus = alien records of Major/Massive species / all
#                             EICAT-matched alien records that year
#   mean_record_eicat       = record-weighted mean EICAT score (MC=0 .. MV=4)
#
# This does NOT guarantee independence from effort - if what gets recorded
# changes composition as effort grows, a ratio will still track effort. That is
# the actual test here.
#
# Writes:
#   report_figures/site_<id>_invasive_impact_composition.csv
#   report_figures/composition_correlations.csv
#   report_figures/fig11_composition_metric.{png,svg,csv}   (new figure)
# Figure 10 and all earlier outputs are left untouched.
# ---------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(impIndicator)

ACCENT_BLUE <- "#1A73E8"
ACCENT_RED <- "#D93025"
ACCENT_GREEN <- "#188038"
NEUTRAL_GREY <- "#70757A"
TITLE_COLOR <- "#202124"
TEXT_COLOR <- "#3C4043"

# Text sizes are set relative to the 9.5 x 10 in canvas: the figure is usually
# embedded at page width, so everything is scaled down on the page. The non-bold
# elements (subtitle, axis text, legend, caption, panel annotations) are sized up
# accordingly - they were unreadable at 7-9 pt.
base_theme <- theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = TITLE_COLOR, hjust = 0,
                              lineheight = 1.15, margin = margin(b = 8)),
    plot.subtitle = element_text(size = 12, color = NEUTRAL_GREY, hjust = 0,
                                 lineheight = 1.2, margin = margin(b = 10)),
    plot.caption = element_text(size = 9.5, color = NEUTRAL_GREY, hjust = 0,
                                lineheight = 1.2, margin = margin(t = 10)),
    axis.title = element_text(face = "bold", size = 12, color = TEXT_COLOR),
    axis.text = element_text(size = 11, color = TEXT_COLOR),
    strip.text = element_text(face = "bold", size = 13, color = TITLE_COLOR, hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 11.5, color = TEXT_COLOR),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.4),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 14, 10, 14)
  )

# Years with fewer alien records than this are kept but flagged, and the
# correlations are reported both with and without them.
LOW_N_FLAG <- 20
EICAT_NUM <- c(MC = 0, MN = 1, MO = 2, MR = 3, MV = 4)
EXPECTED_RAW <- c("786" = 0.798, "67" = 0.990, "879" = 0.942)

gidias_file <- "inst/extdata/GIDIAS/GIDIAS_machine_read.csv"
out_dir <- "report_figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Loading GIDIAS database...\n")
gidias <- read_csv(gidias_file, show_col_types = FALSE)

sites <- list(
  france = list(
    id = 786, name = "La Petite Camargue (France)",
    file = "temp_worker_full_site_786_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/france/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/france/speciesprofile.txt"
  ),
  uk = list(
    id = 67, name = "Severn Estuary (United Kingdom)",
    file = "temp_worker_full_site_67_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/united_kingdom/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/united_kingdom/speciesprofile.txt"
  ),
  botswana = list(
    id = 879, name = "Okavango Delta System (Botswana)",
    file = "temp_worker_full_site_879_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/Botswana/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/Botswana/speciesprofile.txt"
  )
)

# Wilson score interval for a proportion - appropriate for the sparse early
# years where a normal approximation would run off the ends of [0, 1].
wilson_ci <- function(k, n, z = 1.959964) {
  out <- matrix(NA_real_, nrow = length(k), ncol = 2)
  ok <- !is.na(k) & !is.na(n) & n > 0
  p <- k[ok] / n[ok]; nn <- n[ok]
  d <- 1 + z^2 / nn
  centre <- p + z^2 / (2 * nn)
  halfw <- z * sqrt(p * (1 - p) / nn + z^2 / (4 * nn^2))
  out[ok, 1] <- (centre - halfw) / d
  out[ok, 2] <- (centre + halfw) / d
  out
}

corr_row <- function(site_id, site_name, variant, x, effort) {
  keep <- !is.na(x) & !is.na(effort)
  x <- x[keep]; effort <- effort[keep]
  if (length(x) < 3 || sd(x) == 0 || sd(effort) == 0) {
    return(tibble(site_id = site_id, site_name = site_name, variant = variant,
                  pearson_r_with_effort = NA_real_, spearman_rho = NA_real_,
                  p_value = NA_real_, n_years = length(x)))
  }
  ct <- suppressWarnings(cor.test(x, effort, method = "pearson"))
  tibble(site_id = site_id, site_name = site_name, variant = variant,
         pearson_r_with_effort = unname(ct$estimate),
         spearman_rho = suppressWarnings(cor(x, effort, method = "spearman")),
         p_value = ct$p.value, n_years = length(x))
}

fold_var <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2) return(NA_real_)
  max(x) / min(x)
}

site_series <- list()
corr_tbl <- list()
dominance_tbl <- list()
summary_tbl <- list()

for (skey in names(sites)) {
  s <- sites[[skey]]
  cat(sprintf("\n=== Site %d: %s ===\n", s$id, s$name))

  taxon_df <- read.delim(s$griis_taxon, stringsAsFactors = FALSE, sep = "\t", quote = "")
  profile_df <- read.delim(s$griis_profile, stringsAsFactors = FALSE, sep = "\t", quote = "")
  taxon_df <- taxon_df %>% mutate(clean_name = str_extract(scientificName, "^[^ ]+ [^ ]+"))
  alien_registry <- taxon_df %>%
    left_join(profile_df, by = "id") %>%
    select(clean_name, isInvasive) %>%
    distinct(clean_name, .keep_all = TRUE)

  cat("Loading site occurrence cube...\n")
  site_cube <- read_csv(
    s$file,
    col_select = c(mgrscellcode, species, specieskey, year, occurrences),
    show_col_types = FALSE
  )

  # Same alien definition as the previous scripts
  site_mapped <- site_cube %>%
    left_join(alien_registry, by = c("species" = "clean_name")) %>%
    mutate(is_alien = !is.na(isInvasive))
  alien_species_list <- unique(site_mapped$species[site_mapped$is_alien])
  rm(site_mapped)

  eicat_data <- gidias %>%
    filter(Verified.Name.GBIF.Taxon %in% alien_species_list) %>%
    filter(!is.na(magnitude.Nature)) %>%
    mutate(
      scientific_name = Verified.Name.GBIF.Taxon,
      impact_category = case_when(
        magnitude.Nature == 0 ~ "MC", magnitude.Nature == 1 ~ "MN",
        magnitude.Nature == 2 ~ "MO", magnitude.Nature == 3 ~ "MR",
        magnitude.Nature == 4 ~ "MV", TRUE ~ "MC"
      ),
      impact_mechanism = coalesce(mechanism.Nature.clean, "Unknown")
    ) %>%
    select(scientific_name, impact_category, impact_mechanism) %>%
    distinct()

  eicat_species <- unique(eicat_data$scientific_name)
  cat(sprintf("GRIIS alien: %d | matched to GIDIAS: %d\n",
              length(alien_species_list), length(eicat_species)))

  # Max EICAT category per species
  species_max_score <- eicat_data %>%
    mutate(score = unname(EICAT_NUM[impact_category])) %>%
    group_by(scientific_name) %>%
    summarise(max_score = max(score), .groups = "drop")

  cat(sprintf("Species by max EICAT category: %s\n",
              paste(sprintf("%s=%d", names(table(species_max_score$max_score)),
                            as.integer(table(species_max_score$max_score))),
                    collapse = " ")))

  # Raw indicator, kept only as the reference series for the correlation table
  cat("Computing raw EICAT impact indicator (reference)...\n")
  res_imp <- compute_impact_indicator(
    cube = site_cube %>% rename(scientificName = species, cellCode = mgrscellcode,
                                taxonKey = specieskey),
    impact_data = eicat_data, method = "mean_cum", ci_type = "none"
  )
  impact_raw_df <- res_imp$impact %>% select(year, impact_raw = diversity_val)

  annual_effort <- site_cube %>%
    group_by(year) %>%
    summarise(total_occurrences = sum(occurrences), .groups = "drop")

  # --- The composition metric -------------------------------------------
  alien_scored <- site_cube %>%
    filter(species %in% eicat_species) %>%
    left_join(species_max_score, by = c("species" = "scientific_name"))

  composition <- alien_scored %>%
    group_by(year) %>%
    summarise(
      alien_records_eicat = sum(occurrences),
      records_major_plus = sum(occurrences[max_score >= 3]),
      mean_record_eicat = sum(occurrences * max_score) / sum(occurrences),
      n_species_eicat = n_distinct(species),
      n_species_major_plus = n_distinct(species[max_score >= 3]),
      .groups = "drop"
    ) %>%
    mutate(prop_records_major_plus = records_major_plus / alien_records_eicat)

  ci <- wilson_ci(composition$records_major_plus, composition$alien_records_eicat)
  composition$prop_major_plus_lo <- ci[, 1]
  composition$prop_major_plus_hi <- ci[, 2]

  # Is the ratio just one hyper-abundant species? Report the biggest single
  # contributor to the Major/Massive record pool.
  dom <- alien_scored %>%
    filter(max_score >= 3) %>%
    group_by(species) %>%
    summarise(recs = sum(occurrences), .groups = "drop") %>%
    arrange(desc(recs))
  total_mp <- sum(dom$recs)
  dominance_tbl[[skey]] <- tibble(
    site_id = s$id, site_name = s$name,
    n_major_plus_species = nrow(dom),
    top_species = if (nrow(dom) > 0) dom$species[1] else NA_character_,
    top_species_share_of_major_plus = if (nrow(dom) > 0) dom$recs[1] / total_mp else NA_real_,
    top3_share_of_major_plus = if (nrow(dom) > 0) sum(head(dom$recs, 3)) / total_mp else NA_real_
  )

  series <- annual_effort %>%
    left_join(composition, by = "year") %>%
    left_join(impact_raw_df, by = "year") %>%
    mutate(
      alien_records_eicat = coalesce(alien_records_eicat, 0),
      low_n_flag = alien_records_eicat < LOW_N_FLAG,
      data_note = ifelse(alien_records_eicat == 0,
                         "no EICAT-matched alien records; ratio undefined",
                         ifelse(low_n_flag,
                                sprintf("low confidence: %d alien records < %d",
                                        as.integer(alien_records_eicat), LOW_N_FLAG),
                                NA_character_)),
      site_id = s$id, site_name = s$name
    ) %>%
    arrange(year)

  write_csv(
    series %>% select(year, total_occurrences, alien_records_eicat, n_species_eicat,
                      n_species_major_plus, records_major_plus,
                      prop_records_major_plus, prop_major_plus_lo, prop_major_plus_hi,
                      mean_record_eicat, impact_raw, low_n_flag, data_note,
                      site_id, site_name),
    file.path(out_dir, sprintf("site_%d_invasive_impact_composition.csv", s$id))
  )

  well_sampled <- series %>% filter(!low_n_flag)

  corr_tbl[[skey]] <- bind_rows(
    corr_row(s$id, s$name, "raw", series$impact_raw, series$total_occurrences),
    corr_row(s$id, s$name, "prop_records_major_plus",
             series$prop_records_major_plus, series$total_occurrences),
    corr_row(s$id, s$name, "mean_record_eicat",
             series$mean_record_eicat, series$total_occurrences),
    corr_row(s$id, s$name, "prop_records_major_plus_n20",
             well_sampled$prop_records_major_plus, well_sampled$total_occurrences),
    corr_row(s$id, s$name, "mean_record_eicat_n20",
             well_sampled$mean_record_eicat, well_sampled$total_occurrences)
  )

  raw_r <- corr_tbl[[skey]]$pearson_r_with_effort[1]
  cat(sprintf("Sanity: raw r = %.4f (expect %.3f)\n", raw_r, EXPECTED_RAW[[as.character(s$id)]]))
  if (!is.finite(raw_r) || abs(raw_r - EXPECTED_RAW[[as.character(s$id)]]) > 0.01) {
    stop(sprintf("SANITY CHECK FAILED (site %d): raw r = %.4f, expected %.3f.",
                 s$id, raw_r, EXPECTED_RAW[[as.character(s$id)]]))
  }

  summary_tbl[[skey]] <- tibble(
    site_id = s$id, site_name = s$name,
    years_total = nrow(series),
    years_low_n = sum(series$low_n_flag, na.rm = TRUE),
    prop_min = min(series$prop_records_major_plus, na.rm = TRUE),
    prop_median = median(series$prop_records_major_plus, na.rm = TRUE),
    prop_max = max(series$prop_records_major_plus, na.rm = TRUE),
    fold_prop = fold_var(series$prop_records_major_plus),
    fold_raw = fold_var(series$impact_raw),
    mean_eicat_min = min(series$mean_record_eicat, na.rm = TRUE),
    mean_eicat_max = max(series$mean_record_eicat, na.rm = TRUE)
  )

  cat(sprintf("prop_records_major_plus: min %.3f | median %.3f | max %.3f (%d low-n years flagged)\n",
              summary_tbl[[skey]]$prop_min, summary_tbl[[skey]]$prop_median,
              summary_tbl[[skey]]$prop_max, summary_tbl[[skey]]$years_low_n))

  site_series[[skey]] <- series
  rm(site_cube, alien_scored)
  gc(verbose = FALSE)
}

all_corr <- bind_rows(corr_tbl)
write_csv(
  all_corr %>% select(site_id, site_name, variant, pearson_r_with_effort,
                      spearman_rho, p_value, n_years),
  file.path(out_dir, "composition_correlations.csv")
)

all_series <- bind_rows(site_series)
dom_tbl <- bind_rows(dominance_tbl)
sum_tbl <- bind_rows(summary_tbl)

# --- Figure 11 -------------------------------------------------------------
# The proportion is plotted on its own natural 0-1 axis (it is interpretable as
# a percentage), with effort overlaid as a scaled dashed reference line.

site_levels <- vapply(sites, function(x) x$name, character(1), USE.NAMES = FALSE)

scale01 <- function(x) {
  lo <- min(x, na.rm = TRUE); hi <- max(x, na.rm = TRUE)
  if (!is.finite(lo) || !is.finite(hi) || hi == lo) return(rep(NA_real_, length(x)))
  (x - lo) / (hi - lo)
}

plot_df <- all_series %>%
  group_by(site_name) %>%
  mutate(effort_scaled = scale01(total_occurrences)) %>%
  ungroup() %>%
  mutate(site_name = factor(site_name, levels = site_levels))

ann_df <- all_corr %>%
  filter(variant %in% c("prop_records_major_plus", "mean_record_eicat",
                        "prop_records_major_plus_n20")) %>%
  select(site_name, variant, pearson_r_with_effort) %>%
  pivot_wider(names_from = variant, values_from = pearson_r_with_effort) %>%
  mutate(
    label = sprintf(
      "prop Major+ r = %.2f  |  well-sampled years only r = %.2f  |  mean EICAT/record r = %.2f",
      prop_records_major_plus, prop_records_major_plus_n20, mean_record_eicat),
    site_name = factor(site_name, levels = site_levels)
  )

# Title written from the observed result.
r_prop <- all_corr$pearson_r_with_effort[all_corr$variant == "prop_records_major_plus"]
r_prop20 <- all_corr$pearson_r_with_effort[all_corr$variant == "prop_records_major_plus_n20"]
n_sites <- length(unique(all_corr$site_id))
n_coupled <- sum(abs(r_prop20) > 0.5, na.rm = TRUE)

fig_title <- if (n_coupled == 0) {
  paste0("A ratio-based composition metric breaks the link with sampling volume ",
         "at all three sites")
} else if (n_coupled < n_sites) {
  sprintf(paste0("A ratio-based composition metric is largely independent of sampling ",
                 "volume, but still tracks it at %d of %d sites"), n_coupled, n_sites)
} else {
  "Even a ratio-based composition metric tracks sampling volume at every site"
}

fig_subtitle <- paste0(
  "Share of EICAT-assessed alien records belonging to species assessed Major (MR) or ",
  "Massive (MV), with 95% Wilson intervals. Unlike the impact indicator this is a ratio, ",
  "so it has no built-in sample-size dependence and needs no rarefaction. Total occurrences ",
  "(dashed, scaled 0-1) shown as the effort reference; open points mark years with fewer ",
  "than 20 alien records."
)

p11 <- ggplot(plot_df, aes(x = year)) +
  geom_ribbon(aes(ymin = prop_major_plus_lo, ymax = prop_major_plus_hi),
              fill = ACCENT_GREEN, alpha = 0.15) +
  geom_line(aes(y = effort_scaled, color = "Total occurrences (effort, scaled)",
                linetype = "Total occurrences (effort, scaled)"), linewidth = 0.8) +
  geom_line(aes(y = prop_records_major_plus, color = "Proportion of alien records Major/Massive",
                linetype = "Proportion of alien records Major/Massive"), linewidth = 0.9) +
  geom_point(aes(y = prop_records_major_plus, shape = low_n_flag),
             color = ACCENT_GREEN, size = 1.5) +
  geom_text(data = ann_df, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.02, vjust = 1.4, size = 3.9, color = TEXT_COLOR,
            family = "sans", fontface = "italic", inherit.aes = FALSE) +
  facet_wrap(~site_name, ncol = 1) +
  # One merged legend: the linetype scale is suppressed and its information is
  # carried into the colour keys via override.aes.
  scale_color_manual(
    values = c("Proportion of alien records Major/Massive" = ACCENT_GREEN,
               "Total occurrences (effort, scaled)" = NEUTRAL_GREY),
    breaks = c("Proportion of alien records Major/Massive",
               "Total occurrences (effort, scaled)")) +
  scale_linetype_manual(
    values = c("Proportion of alien records Major/Massive" = "solid",
               "Total occurrences (effort, scaled)" = "dashed"),
    guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), guide = "none") +
  # Extra headroom so the larger annotation clears the data.
  scale_y_continuous(limits = c(0, 1.22), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  labs(
    # Wrap widths are tied to the font sizes above - raising the sizes without
    # narrowing these would run the text off the canvas.
    title = str_wrap(fig_title, 70),
    subtitle = str_wrap(fig_subtitle, 88),
    x = "Year", y = "Proportion of alien records Major/Massive",
    caption = str_wrap(paste0(
      "Source: GBIF occurrence cubes for three Ramsar sites; GRIIS national alien checklists; ",
      "GIDIAS EICAT assessments. Annotations give the Pearson correlation of each series with ",
      "annual total occurrences. Effort line is scaled 0-1 within site and shares the axis for ",
      "shape comparison only."
    ), 125)
  ) +
  base_theme +
  guides(color = guide_legend(
    nrow = 1,
    override.aes = list(linetype = c("solid", "dashed"), shape = NA)))

# Taller canvas at the same width: the header/legend/caption block grew with the
# larger type, so this gives the panels their height back. Embedding scales by
# width, so the on-page text size is unaffected.
ggsave(file.path(out_dir, "fig11_composition_metric.png"),
       plot = p11, width = 9.5, height = 12.0, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig11_composition_metric.svg"),
       plot = p11, width = 9.5, height = 12.0, bg = "white")

write_csv(
  plot_df %>% select(site_id, site_name, year, total_occurrences, effort_scaled,
                     alien_records_eicat, prop_records_major_plus,
                     prop_major_plus_lo, prop_major_plus_hi, mean_record_eicat,
                     low_n_flag),
  file.path(out_dir, "fig11_composition_metric.csv")
)

# --- Console summary -------------------------------------------------------

cat("\n\n=============================================================\n")
cat("COMPOSITION METRIC - CORRELATION WITH TOTAL OCCURRENCES\n")
cat("=============================================================\n")
print(as.data.frame(all_corr %>%
  mutate(pearson_r_with_effort = round(pearson_r_with_effort, 3),
         spearman_rho = round(spearman_rho, 3),
         p_value = signif(p_value, 3))), row.names = FALSE)

cat("\n--- Series summary ---\n")
print(as.data.frame(sum_tbl %>%
  mutate(across(c(prop_min, prop_median, prop_max, mean_eicat_min, mean_eicat_max),
                ~round(.x, 3)),
         across(starts_with("fold_"), ~round(.x, 1)))), row.names = FALSE)

cat("\n--- Is the ratio driven by one species? ---\n")
print(as.data.frame(dom_tbl %>%
  mutate(across(ends_with("_share_of_major_plus"), ~round(.x, 3)))), row.names = FALSE)

cat("\n--- Files written ---\n")
cat("  report_figures/site_<id>_invasive_impact_composition.csv (3 files)\n")
cat("  report_figures/composition_correlations.csv\n")
cat("  report_figures/fig11_composition_metric.png / .svg / .csv\n")
cat("\nFigure title selected from the observed result:\n  ", fig_title, "\n", sep = "")
