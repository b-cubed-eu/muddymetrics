# ---------------------------------------------------------------------------
# Effort-standardised invasive-impact indicators (Part 3 extension)
#
# Question: once sampling effort is controlled, does an interpretable invasive
# impact signal remain, or does it disappear entirely?
#
# Variant A - rarefaction to a constant annual record count.
# Variant B - occupancy-based (presence per 100 m MGRS cell).
#
# Reuses the data-loading logic of scripts/generate_all_sites_invasive.R.
# Writes:
#   report_figures/site_<id>_invasive_impact_standardised.csv
#   report_figures/effort_standardisation_correlations.csv
#   report_figures/fig10a_size_rarefied_impact.{png,svg,csv}   (with occupancy)
#   report_figures/fig10b_size_rarefied_only.{png,svg,csv}     (rarefaction only)
#
# NOTE ON THE FIGURE NAME: this figure was originally written as
# fig10_effort_standardised_impact.*. That name now belongs to the
# coverage-based analysis (scripts/generate_coverage_standardised_invasive.R),
# so this one writes to fig10a_* instead. Both figures therefore coexist and
# re-running either script no longer overwrites the other.
# ---------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(tidyr)
library(impIndicator)

# Sizing & Theme Constants (matching the existing report figures)
ACCENT_BLUE <- "#1A73E8"
ACCENT_RED <- "#D93025"
ACCENT_GREEN <- "#188038"
NEUTRAL_GREY <- "#70757A"
TITLE_COLOR <- "#202124"
TEXT_COLOR <- "#3C4043"

base_theme <- theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 13, color = TITLE_COLOR, hjust = 0,
                              lineheight = 1.15, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 9, color = NEUTRAL_GREY, hjust = 0),
    plot.caption = element_text(size = 7, color = NEUTRAL_GREY, hjust = 0),
    axis.title = element_text(face = "bold", size = 9, color = TEXT_COLOR),
    axis.text = element_text(size = 8, color = TEXT_COLOR),
    strip.text = element_text(face = "bold", size = 9.5, color = TITLE_COLOR, hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 8, color = TEXT_COLOR),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E0E0E0", linewidth = 0.4),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(12, 14, 10, 14)
  )

# --- Analysis constants ----------------------------------------------------
N_MIN_FLOOR <- 20    # years with fewer EICAT-matched alien records are excluded
B_DRAWS <- 100       # rarefaction draws
EICAT_NUM <- c(MC = 0, MN = 1, MO = 2, MR = 3, MV = 4)

gidias_file <- "inst/extdata/GIDIAS/GIDIAS_machine_read.csv"
out_dir <- "report_figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("Loading GIDIAS database...\n")
gidias <- read_csv(gidias_file, show_col_types = FALSE)

sites <- list(
  france = list(
    id = 786,
    name = "La Petite Camargue (France)",
    file = "temp_worker_full_site_786_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/france/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/france/speciesprofile.txt",
    expected_raw_r = 0.798
  ),
  uk = list(
    id = 67,
    name = "Severn Estuary (United Kingdom)",
    file = "temp_worker_full_site_67_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/united_kingdom/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/united_kingdom/speciesprofile.txt",
    expected_raw_r = 0.990
  ),
  botswana = list(
    id = 879,
    name = "Okavango Delta System (Botswana)",
    file = "temp_worker_full_site_879_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/Botswana/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/Botswana/speciesprofile.txt",
    expected_raw_r = 0.942
  )
)

# --- Helpers ---------------------------------------------------------------

# Correlate one series against sampling effort; returns a one-row tibble.
corr_row <- function(site_id, site_name, variant, x, effort) {
  keep <- !is.na(x) & !is.na(effort)
  x <- x[keep]; effort <- effort[keep]
  if (length(x) < 3 || sd(x) == 0 || sd(effort) == 0) {
    return(tibble(
      site_id = site_id, site_name = site_name, variant = variant,
      pearson_r_with_effort = NA_real_, spearman_rho = NA_real_,
      p_value = NA_real_, n_years = length(x)
    ))
  }
  ct <- suppressWarnings(cor.test(x, effort, method = "pearson"))
  rho <- suppressWarnings(cor(x, effort, method = "spearman"))
  tibble(
    site_id = site_id, site_name = site_name, variant = variant,
    pearson_r_with_effort = unname(ct$estimate), spearman_rho = rho,
    p_value = ct$p.value, n_years = length(x)
  )
}

# Scale a vector to 0-1 using a supplied reference range (for the CI band we
# reuse the range of the mean series so the band stays on the same scale).
scale01 <- function(x, lo = min(x, na.rm = TRUE), hi = max(x, na.rm = TRUE)) {
  if (!is.finite(lo) || !is.finite(hi) || hi == lo) return(rep(NA_real_, length(x)))
  (x - lo) / (hi - lo)
}

site_series <- list()
corr_tbl <- list()
rarefaction_log <- list()

# --- Main loop -------------------------------------------------------------

for (skey in names(sites)) {
  s <- sites[[skey]]
  cat(sprintf("\n=== Site %d: %s ===\n", s$id, s$name))

  # 1. GRIIS alien registry (same logic as generate_all_sites_invasive.R)
  taxon_df <- read.delim(s$griis_taxon, stringsAsFactors = FALSE, sep = "\t", quote = "")
  profile_df <- read.delim(s$griis_profile, stringsAsFactors = FALSE, sep = "\t", quote = "")

  taxon_df <- taxon_df %>%
    mutate(clean_name = str_extract(scientificName, "^[^ ]+ [^ ]+"))

  alien_registry <- taxon_df %>%
    left_join(profile_df, by = "id") %>%
    select(clean_name, isInvasive) %>%
    distinct(clean_name, .keep_all = TRUE)

  # 2. Occurrence cube (only the columns we need, these files are ~150 MB)
  cat("Loading site occurrence cube...\n")
  site_cube <- read_csv(
    s$file,
    col_select = c(mgrscellcode, species, specieskey, year, occurrences),
    show_col_types = FALSE
  )

  site_mapped <- site_cube %>%
    left_join(alien_registry, by = c("species" = "clean_name")) %>%
    mutate(is_alien = !is.na(isInvasive))

  alien_species_list <- unique(site_mapped$species[site_mapped$is_alien])
  cat(sprintf("Site species: %d | GRIIS alien: %d\n",
              length(unique(site_cube$species)), length(alien_species_list)))

  # 3. GIDIAS EICAT match (identical construction to the Part 3 pilot)
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

  eicat_species <- unique(eicat_data$scientific_name)
  cat(sprintf("Matched to GIDIAS: %d alien species\n", length(eicat_species)))

  # Max EICAT category per species (used by Variant B; unmatched = not assessed)
  species_max_score <- eicat_data %>%
    mutate(score = unname(EICAT_NUM[impact_category])) %>%
    group_by(scientific_name) %>%
    summarise(max_score = max(score), .groups = "drop")

  # 4. Raw indicator on the full cube -------------------------------------
  site_cube_imp <- site_cube %>%
    rename(scientificName = species, cellCode = mgrscellcode, taxonKey = specieskey)

  cat("Computing raw EICAT impact indicator...\n")
  res_imp <- compute_impact_indicator(
    cube = site_cube_imp,
    impact_data = eicat_data,
    method = "mean_cum",
    ci_type = "none"
  )
  impact_raw_df <- res_imp$impact %>%
    select(year, impact_raw = diversity_val)

  # Annual effort and EICAT-matched alien record totals
  annual_effort <- site_cube %>%
    group_by(year) %>%
    summarise(total_occurrences = sum(occurrences), .groups = "drop")

  alien_rec <- site_cube %>%
    filter(species %in% eicat_species) %>%
    transmute(
      cellCode = mgrscellcode, scientificName = species,
      taxonKey = specieskey, year, occurrences
    )

  annual_alien <- alien_rec %>%
    group_by(year) %>%
    summarise(alien_records_eicat = sum(occurrences), .groups = "drop")

  # 5. Variant A - rarefaction --------------------------------------------
  excluded <- annual_alien %>% filter(alien_records_eicat < N_MIN_FLOOR)
  retained <- annual_alien %>% filter(alien_records_eicat >= N_MIN_FLOOR)
  # Years present in the cube but with zero EICAT-matched alien records
  zero_years <- setdiff(annual_effort$year, annual_alien$year)

  n_rare <- if (nrow(retained) > 0) min(retained$alien_records_eicat) else NA_integer_

  cat(sprintf("Rarefaction: n = %s records/year | retained years = %d | excluded years = %d\n",
              ifelse(is.na(n_rare), "NA", format(n_rare)),
              nrow(retained), nrow(excluded) + length(zero_years)))
  if (nrow(excluded) > 0) {
    cat("  Excluded (below n_min = 20): ",
        paste(sprintf("%d (%d rec)", excluded$year, excluded$alien_records_eicat),
              collapse = ", "), "\n", sep = "")
  }
  if (length(zero_years) > 0) {
    cat("  Excluded (no EICAT-matched alien records): ",
        paste(sort(zero_years), collapse = ", "), "\n", sep = "")
  }

  rarefied_summary <- tibble(
    year = integer(0), impact_rarefied_mean = numeric(0),
    impact_rarefied_lo = numeric(0), impact_rarefied_hi = numeric(0)
  )

  if (!is.na(n_rare) && nrow(retained) >= 3) {
    rare_years <- sort(retained$year)

    # Pre-expand each retained year to record level once: a row with k
    # occurrences contributes k draws, so it is k times more likely to be hit.
    year_rows <- lapply(rare_years, function(y) alien_rec[alien_rec$year == y, ])
    names(year_rows) <- as.character(rare_years)
    year_idx <- lapply(year_rows, function(d) rep.int(seq_len(nrow(d)), d$occurrences))

    set.seed(20260729L + s$id)
    cat(sprintf("Running %d rarefaction draws", B_DRAWS))
    draw_mat <- matrix(NA_real_, nrow = length(rare_years), ncol = B_DRAWS,
                       dimnames = list(as.character(rare_years), NULL))

    for (b in seq_len(B_DRAWS)) {
      drawn <- lapply(seq_along(rare_years), function(i) {
        d <- year_rows[[i]]
        idx <- year_idx[[i]]
        hit <- tabulate(sample(idx, n_rare, replace = FALSE), nbins = nrow(d))
        keep <- hit > 0
        d_sub <- d[keep, c("cellCode", "scientificName", "taxonKey", "year")]
        d_sub$occurrences <- hit[keep]
        d_sub
      })
      cube_b <- bind_rows(drawn)

      res_b <- suppressMessages(compute_impact_indicator(
        cube = cube_b, impact_data = eicat_data,
        method = "mean_cum", ci_type = "none"
      ))
      vals <- res_b$impact
      draw_mat[as.character(vals$year), b] <- vals$diversity_val
      if (b %% 10 == 0) cat(".")
    }
    cat(" done\n")

    # Years absent from a draw contribute a genuine zero, not a missing value.
    draw_mat[is.na(draw_mat)] <- 0

    rarefied_summary <- tibble(
      year = as.integer(rownames(draw_mat)),
      impact_rarefied_mean = rowMeans(draw_mat),
      impact_rarefied_lo = apply(draw_mat, 1, quantile, probs = 0.025, names = FALSE),
      impact_rarefied_hi = apply(draw_mat, 1, quantile, probs = 0.975, names = FALSE)
    )
  } else {
    cat("  Too few retained years for rarefaction - Variant A skipped.\n")
  }

  # 6. Variant B - occupancy ----------------------------------------------
  occupied_cells <- site_cube %>% distinct(year, mgrscellcode)

  cell_alien <- site_cube %>%
    filter(species %in% eicat_species) %>%
    left_join(species_max_score, by = c("species" = "scientific_name")) %>%
    group_by(year, mgrscellcode) %>%
    summarise(cell_max_score = max(max_score), .groups = "drop")

  occupancy <- occupied_cells %>%
    left_join(cell_alien, by = c("year", "mgrscellcode")) %>%
    mutate(cell_max_score = coalesce(cell_max_score, 0)) %>%
    group_by(year) %>%
    summarise(
      n_cells_occupied = n(),
      prop_cells_major_plus = mean(cell_max_score >= 3),
      mean_cell_eicat = mean(cell_max_score),
      .groups = "drop"
    )

  # 7. Assemble the per-site series ---------------------------------------
  excl_reason <- bind_rows(
    excluded %>% transmute(
      year,
      year_excluded_reason = sprintf("excluded from rarefaction: %d EICAT-matched alien records < n_min = %d",
                                     alien_records_eicat, N_MIN_FLOOR)
    ),
    tibble(year = zero_years,
           year_excluded_reason = "excluded from rarefaction: no EICAT-matched alien records")
  )

  series <- annual_effort %>%
    left_join(annual_alien, by = "year") %>%
    mutate(alien_records_eicat = coalesce(alien_records_eicat, 0L)) %>%
    left_join(impact_raw_df, by = "year") %>%
    left_join(rarefied_summary, by = "year") %>%
    left_join(excl_reason, by = "year") %>%
    left_join(occupancy, by = "year") %>%
    mutate(
      rarefaction_n = ifelse(is.na(impact_rarefied_mean), NA_integer_, as.integer(n_rare)),
      year_excluded_reason = coalesce(year_excluded_reason, NA_character_),
      site_id = s$id,
      site_name = s$name
    ) %>%
    arrange(year) %>%
    select(year, total_occurrences, alien_records_eicat, impact_raw,
           impact_rarefied_mean, impact_rarefied_lo, impact_rarefied_hi,
           rarefaction_n, year_excluded_reason,
           n_cells_occupied, prop_cells_major_plus, mean_cell_eicat,
           site_id, site_name)

  write_csv(series, file.path(out_dir, sprintf("site_%d_invasive_impact_standardised.csv", s$id)))

  # 8. Correlations vs effort ---------------------------------------------
  corr_tbl[[skey]] <- bind_rows(
    corr_row(s$id, s$name, "raw", series$impact_raw, series$total_occurrences),
    corr_row(s$id, s$name, "rarefied", series$impact_rarefied_mean, series$total_occurrences),
    corr_row(s$id, s$name, "occupancy_prop_major", series$prop_cells_major_plus, series$total_occurrences),
    corr_row(s$id, s$name, "occupancy_mean_eicat", series$mean_cell_eicat, series$total_occurrences)
  )

  raw_r <- corr_tbl[[skey]]$pearson_r_with_effort[1]
  cat(sprintf("Raw r vs effort = %.4f (pilot reported %.3f)\n", raw_r, s$expected_raw_r))
  if (!is.finite(raw_r) || abs(raw_r - s$expected_raw_r) > 0.01) {
    stop(sprintf(
      "SANITY CHECK FAILED for site %d: raw r = %.4f but the Part 3 pilot reported %.3f. Stopping.",
      s$id, raw_r, s$expected_raw_r
    ))
  }

  rarefaction_log[[skey]] <- tibble(
    site_id = s$id, site_name = s$name,
    rarefaction_n = as.integer(n_rare),
    years_retained = nrow(retained),
    years_excluded = nrow(excluded) + length(zero_years),
    excluded_years = paste(sort(c(excluded$year, zero_years)), collapse = "; ")
  )

  site_series[[skey]] <- series
  rm(site_cube, site_mapped, site_cube_imp, alien_rec)
  gc(verbose = FALSE)
}

all_corr <- bind_rows(corr_tbl)
write_csv(
  all_corr %>% select(site_id, site_name, variant, pearson_r_with_effort,
                      spearman_rho, p_value, n_years),
  file.path(out_dir, "effort_standardisation_correlations.csv")
)

all_series <- bind_rows(site_series)
rare_log <- bind_rows(rarefaction_log)

# --- Figure 10 -------------------------------------------------------------

site_levels <- vapply(sites, function(x) x$name, character(1), USE.NAMES = FALSE)

series_long <- all_series %>%
  group_by(site_name) %>%
  mutate(
    `Raw impact (EICAT mean_cum)` = scale01(impact_raw),
    `Rarefied impact (constant effort)` = scale01(impact_rarefied_mean),
    `Occupancy (prop. cells Major+)` = scale01(prop_cells_major_plus),
    `Total occurrences (effort)` = scale01(total_occurrences)
  ) %>%
  ungroup() %>%
  select(site_name, year,
         `Raw impact (EICAT mean_cum)`, `Rarefied impact (constant effort)`,
         `Occupancy (prop. cells Major+)`, `Total occurrences (effort)`) %>%
  pivot_longer(-c(site_name, year), names_to = "series", values_to = "value_scaled") %>%
  # NA rows are kept deliberately: years excluded from rarefaction must appear
  # as gaps in the line, not be bridged over.
  mutate(site_name = factor(site_name, levels = site_levels))

ribbon_df <- all_series %>%
  group_by(site_name) %>%
  mutate(
    lo_scaled = scale01(impact_rarefied_lo,
                        min(impact_rarefied_mean, na.rm = TRUE),
                        max(impact_rarefied_mean, na.rm = TRUE)),
    hi_scaled = scale01(impact_rarefied_hi,
                        min(impact_rarefied_mean, na.rm = TRUE),
                        max(impact_rarefied_mean, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  filter(!is.na(impact_rarefied_mean)) %>%
  select(site_name, year, lo_scaled, hi_scaled) %>%
  mutate(site_name = factor(site_name, levels = site_levels))

ann_df <- all_corr %>%
  select(site_name, variant, pearson_r_with_effort) %>%
  pivot_wider(names_from = variant, values_from = pearson_r_with_effort) %>%
  mutate(
    label = sprintf("raw r = %.2f  |  rarefied r = %.2f  |  occupancy r = %.2f",
                    raw, rarefied, occupancy_prop_major),
    site_name = factor(site_name, levels = site_levels)
  )

series_levels <- c("Raw impact (EICAT mean_cum)", "Rarefied impact (constant effort)",
                   "Occupancy (prop. cells Major+)", "Total occurrences (effort)")
series_colors <- setNames(c(ACCENT_RED, ACCENT_BLUE, ACCENT_GREEN, NEUTRAL_GREY), series_levels)
series_ltys <- setNames(c("solid", "solid", "solid", "dashed"), series_levels)
series_shapes <- setNames(c(16, 17, 15, NA), series_levels)

series_long$series <- factor(series_long$series, levels = series_levels)

# Title states the finding actually observed. "Still effort-coupled" means the
# strong positive correlation seen in the raw series survives standardisation.
r_of <- function(v) all_corr$pearson_r_with_effort[all_corr$variant == v]
n_sites <- length(unique(all_corr$site_id))
rar_still <- sum(r_of("rarefied") > 0.5, na.rm = TRUE)
occ_still <- sum(r_of("occupancy_prop_major") > 0.5, na.rm = TRUE)

fig_title <- if (rar_still == 0 && occ_still == 0) {
  "Standardising for effort decouples invasive impact from sampling volume"
} else if (rar_still == 0 && occ_still > 0) {
  sprintf(paste0("Rarefaction breaks the link between invasive impact and sampling volume, ",
                 "but occupancy still tracks effort at %d of %d sites"), occ_still, n_sites)
} else if (occ_still == 0) {
  sprintf(paste0("Occupancy standardisation decouples invasive impact from sampling volume, ",
                 "but the rarefied series still tracks effort at %d of %d sites"), rar_still, n_sites)
} else {
  "Invasive-impact indicators stay tied to sampling volume even after effort standardisation"
}

fig_subtitle <- sprintf(
  paste0("Each series scaled 0-1 within its own site and metric so shapes are comparable. ",
         "Rarefaction n per site: %s (B = %d draws, shaded 95%% band); occupancy is the ",
         "proportion of occupied 100 m cells holding a Major/Massive alien."),
  paste(sprintf("%s = %d", rare_log$site_id, rare_log$rarefaction_n), collapse = ", "),
  B_DRAWS
)

p10 <- ggplot() +
  geom_ribbon(data = ribbon_df,
              aes(x = year, ymin = lo_scaled, ymax = hi_scaled),
              fill = ACCENT_BLUE, alpha = 0.15) +
  geom_line(data = series_long,
            aes(x = year, y = value_scaled, color = series, linetype = series),
            linewidth = 0.8) +
  geom_point(data = series_long %>% filter(series != "Total occurrences (effort)"),
             aes(x = year, y = value_scaled, color = series, shape = series),
             size = 1.3) +
  geom_text(data = ann_df, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.03, vjust = 1.5, size = 2.9, color = TEXT_COLOR,
            family = "sans", fontface = "italic") +
  facet_wrap(~site_name, ncol = 1) +
  scale_color_manual(values = series_colors, breaks = series_levels) +
  scale_linetype_manual(values = series_ltys, breaks = series_levels) +
  scale_shape_manual(values = series_shapes, breaks = series_levels, guide = "none") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  # Clip rather than drop, so the wide rarefaction band stays continuous.
  coord_cartesian(ylim = c(-0.62, 1.45)) +
  labs(
    title = str_wrap(fig_title, 88),
    subtitle = str_wrap(fig_subtitle, 118),
    x = "Year", y = "Scaled value (0-1 within series)",
    caption = paste0(
      "Source: GBIF occurrence cubes for three Ramsar sites; GRIIS national alien checklists; ",
      "GIDIAS EICAT assessments. Indicator: impIndicator::compute_impact_indicator(method = \"mean_cum\").\n",
      "Annotations give the Pearson correlation of each series with annual total occurrences."
    )
  ) +
  base_theme +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

ggsave(file.path(out_dir, "fig10a_size_rarefied_impact.png"),
       plot = p10, width = 9.5, height = 10.0, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig10a_size_rarefied_impact.svg"),
       plot = p10, width = 9.5, height = 10.0, bg = "white")

fig_csv <- series_long %>%
  mutate(series_label = as.character(series), .keep = "unused") %>%
  left_join(
    ribbon_df %>% mutate(series_label = "Rarefied impact (constant effort)"),
    by = c("site_name", "year", "series_label")
  ) %>%
  arrange(site_name, series_label, year)
write_csv(fig_csv, file.path(out_dir, "fig10a_size_rarefied_impact.csv"))

# --- Figure 10b: rarefaction only (no occupancy) ---------------------------
# Same data and styling as 10a with the occupancy series dropped, for reporting
# the size-based rarefaction result on its own.

fold_var <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2) return(NA_real_)
  max(x) / min(x)
}

fold_b <- all_series %>%
  group_by(site_id, site_name) %>%
  summarise(fold_raw = fold_var(impact_raw),
            fold_rarefied = fold_var(impact_rarefied_mean), .groups = "drop")

levels_b <- c("Raw impact (EICAT mean_cum)", "Rarefied impact (constant effort)",
              "Total occurrences (effort)")
series_b <- series_long %>% filter(series %in% levels_b) %>%
  mutate(series = factor(as.character(series), levels = levels_b))

ann_b <- all_corr %>%
  filter(variant %in% c("raw", "rarefied")) %>%
  select(site_name, variant, pearson_r_with_effort) %>%
  pivot_wider(names_from = variant, values_from = pearson_r_with_effort) %>%
  left_join(fold_b %>% select(site_name, fold_raw, fold_rarefied), by = "site_name") %>%
  mutate(
    label = sprintf("raw r = %.2f (%.0fx variation)  |  rarefied r = %.2f (%.1fx variation)",
                    raw, fold_raw, rarefied, fold_rarefied),
    site_name = factor(site_name, levels = site_levels)
  )

# Title written from the observed result: does rarefaction remove the positive
# effort coupling, and what is left of the indicator's variation once it does?
rar_r <- all_corr$pearson_r_with_effort[all_corr$variant == "rarefied"]
rar_pos_gone <- all(rar_r < 0.3, na.rm = TRUE)
rar_flat <- all(fold_b$fold_rarefied < 0.1 * fold_b$fold_raw, na.rm = TRUE)

title_b <- if (rar_pos_gone && rar_flat) {
  paste0("Holding record count constant removes the effort correlation - and almost all ",
         "of the indicator's variation with it")
} else if (rar_pos_gone) {
  "Rarefaction breaks the link between invasive impact and sampling volume"
} else {
  sprintf("Rarefied invasive impact still tracks sampling volume at %d of %d sites",
          sum(rar_r > 0.5, na.rm = TRUE), length(unique(all_corr$site_id)))
}

subtitle_b <- sprintf(
  paste0("Impact rarefied to a constant %s records per year (B = %d draws, shaded 95%% ",
         "band); series scaled 0-1 within site and metric. Fold-variation (max/min) is ",
         "given unscaled in each panel: the raw indicator moves %.0f-%.0fx across the ",
         "series, the rarefied one %.1f-%.1fx. Gaps are years excluded for holding fewer ",
         "than %d EICAT-matched alien records."),
  paste(sprintf("%s (site %s)", rare_log$rarefaction_n, rare_log$site_id), collapse = ", "),
  B_DRAWS, min(fold_b$fold_raw), max(fold_b$fold_raw),
  min(fold_b$fold_rarefied), max(fold_b$fold_rarefied), N_MIN_FLOOR
)

p10b <- ggplot() +
  geom_ribbon(data = ribbon_df, aes(x = year, ymin = lo_scaled, ymax = hi_scaled),
              fill = ACCENT_BLUE, alpha = 0.15) +
  geom_line(data = series_b,
            aes(x = year, y = value_scaled, color = series, linetype = series),
            linewidth = 0.8) +
  geom_point(data = series_b %>% filter(series != "Total occurrences (effort)"),
             aes(x = year, y = value_scaled, color = series, shape = series), size = 1.3) +
  geom_text(data = ann_b, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.03, vjust = 1.5, size = 2.8, color = TEXT_COLOR,
            family = "sans", fontface = "italic") +
  facet_wrap(~site_name, ncol = 1) +
  scale_color_manual(values = series_colors[levels_b], breaks = levels_b) +
  scale_linetype_manual(values = series_ltys[levels_b], breaks = levels_b) +
  scale_shape_manual(values = series_shapes[levels_b], breaks = levels_b, guide = "none") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  coord_cartesian(ylim = c(-0.62, 1.45)) +
  labs(
    title = str_wrap(title_b, 88),
    subtitle = str_wrap(subtitle_b, 118),
    x = "Year", y = "Scaled value (0-1 within series)",
    caption = paste0(
      "Source: GBIF occurrence cubes for three Ramsar sites; GRIIS national alien checklists; ",
      "GIDIAS EICAT assessments. Indicator: impIndicator::compute_impact_indicator(method = \"mean_cum\").\n",
      "Annotations give the Pearson correlation of each series with annual total occurrences."
    )
  ) +
  base_theme +
  guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1))

ggsave(file.path(out_dir, "fig10b_size_rarefied_only.png"),
       plot = p10b, width = 9.5, height = 10.0, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig10b_size_rarefied_only.svg"),
       plot = p10b, width = 9.5, height = 10.0, bg = "white")

write_csv(
  series_b %>%
    mutate(series_label = as.character(series), .keep = "unused") %>%
    left_join(ribbon_df %>% mutate(series_label = "Rarefied impact (constant effort)"),
              by = c("site_name", "year", "series_label")) %>%
    arrange(site_name, series_label, year),
  file.path(out_dir, "fig10b_size_rarefied_only.csv")
)

# --- Console summary -------------------------------------------------------

cat("\n\n=============================================================\n")
cat("EFFORT STANDARDISATION - CORRELATION WITH TOTAL OCCURRENCES\n")
cat("=============================================================\n")
print(as.data.frame(
  all_corr %>%
    mutate(
      pearson_r_with_effort = round(pearson_r_with_effort, 3),
      spearman_rho = round(spearman_rho, 3),
      p_value = signif(p_value, 3)
    ) %>%
    select(site_id, site_name, variant, pearson_r_with_effort, spearman_rho, p_value, n_years)
), row.names = FALSE)

cat("\n--- Rarefaction settings actually used ---\n")
print(as.data.frame(rare_log), row.names = FALSE)

cat("\n--- Figures / tables written ---\n")
cat(" ", file.path(out_dir, "site_<id>_invasive_impact_standardised.csv"), "(3 files)\n")
cat(" ", file.path(out_dir, "effort_standardisation_correlations.csv"), "\n")
cat(" ", file.path(out_dir, "fig10a_size_rarefied_impact.png / .svg / .csv"), "\n")
cat(" ", file.path(out_dir, "fig10b_size_rarefied_only.png / .svg / .csv"), "(no occupancy)\n")
cat("\nFigure title selected from the observed result:\n  ", fig_title, "\n", sep = "")
