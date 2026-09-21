# ---------------------------------------------------------------------------
# Coverage-based rarefaction for the EICAT impact indicator (Part 3 extension)
#
# Supersedes the size-based rarefaction in
# scripts/generate_effort_standardised_invasive.R. Chao & Jost (2012): equal
# sample SIZE does not imply equal sample COMPLETENESS when assemblages differ
# in evenness, so we standardise to a common estimated coverage C* instead of a
# common record count.
#
# Both rarefactions are recomputed here so the v2 correlation file is
# self-contained and the size-based result acts as an end-to-end sanity check.
# The occupancy variants are deliberately NOT recomputed.
#
# Writes:
#   report_figures/site_<id>_invasive_impact_coverage.csv
#   report_figures/effort_standardisation_correlations_v2.csv
#   report_figures/fig10_effort_standardised_impact.{png,svg,csv}  (replaced)
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

# --- Analysis constants ----------------------------------------------------
N_MIN_FLOOR <- 20     # size-based rarefaction: minimum records for a year to count
COVERAGE_FLOOR <- 0.5 # coverage-based: years below this are too incomplete
B_DRAWS <- 100

# Expected r values from the size-based run (sanity checks)
EXPECTED_RAW <- c("786" = 0.798, "67" = 0.990, "879" = 0.942)
EXPECTED_SIZE <- c("786" = -0.310, "67" = -0.620, "879" = -0.421)

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
    griis_profile = "inst/extdata/griis_checklists/france/speciesprofile.txt"
  ),
  uk = list(
    id = 67,
    name = "Severn Estuary (United Kingdom)",
    file = "temp_worker_full_site_67_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/united_kingdom/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/united_kingdom/speciesprofile.txt"
  ),
  botswana = list(
    id = 879,
    name = "Okavango Delta System (Botswana)",
    file = "temp_worker_full_site_879_full.csv",
    griis_taxon = "inst/extdata/griis_checklists/Botswana/taxon.txt",
    griis_profile = "inst/extdata/griis_checklists/Botswana/speciesprofile.txt"
  )
)

# --- Coverage machinery ----------------------------------------------------

# Chao's abundance-based sample-coverage estimator for the observed sample:
#   C_hat = 1 - (f1/n) * [ (n-1)f1 / ((n-1)f1 + 2f2) ]
# x is a vector of per-species record counts for one year.
coverage_hat <- function(x) {
  x <- x[x > 0]
  n <- sum(x)
  if (n < 2) return(NA_real_)          # undefined with fewer than 2 records
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  if (f1 == 0) return(1)               # no singletons -> coverage 1, and avoids 0/0
  a <- (n - 1) * f1
  1 - (f1 / n) * (a / (a + 2 * f2))
}

# Expected coverage of a size-m subsample drawn without replacement (Chao &
# Jost 2012 interpolation). At m = n the interpolation degenerates to 1, so the
# observed-sample estimator above is used for the endpoint instead.
coverage_at_m <- function(x, m, c_obs) {
  x <- x[x > 0]
  n <- sum(x)
  if (m >= n) return(c_obs)
  if (m < 1) return(0)
  xx <- x[(n - x) >= m]
  if (length(xx) == 0) return(1)
  1 - sum(xx / n * exp(lgamma(n - xx + 1) - lgamma(n - xx - m + 1) -
                         lgamma(n) + lgamma(n - m)))
}

# Smallest m in 1..n whose expected coverage reaches c_target. Coverage is
# monotone increasing in m, so a binary search is valid. A solution always
# exists because c_target is the minimum observed coverage across retained
# years, so m = n always satisfies it.
solve_m <- function(x, c_target, c_obs) {
  n <- sum(x[x > 0])
  if (coverage_at_m(x, 1, c_obs) >= c_target) return(1L)
  lo <- 1L; hi <- as.integer(n)
  while (lo + 1L < hi) {
    mid <- (lo + hi) %/% 2L
    if (coverage_at_m(x, mid, c_obs) >= c_target) hi <- mid else lo <- mid
  }
  hi
}

# --- Rarefaction driver ----------------------------------------------------
# year_rows / year_idx: per-year cube rows and their record-level expansion (a
# row with k occurrences appears k times, so it is k times more likely to be
# drawn). m_by_year: named vector of draw sizes, names are years, in the order
# the draws should be made.
run_rarefaction <- function(year_rows, year_idx, m_by_year, eicat_data, seed, tag) {
  yrs <- names(m_by_year)
  set.seed(seed)
  cat(sprintf("Running %d %s draws", B_DRAWS, tag))
  mat <- matrix(NA_real_, nrow = length(yrs), ncol = B_DRAWS,
                dimnames = list(yrs, NULL))

  for (b in seq_len(B_DRAWS)) {
    drawn <- lapply(yrs, function(y) {
      d <- year_rows[[y]]
      idx <- year_idx[[y]]
      hit <- tabulate(sample(idx, m_by_year[[y]], replace = FALSE), nbins = nrow(d))
      keep <- hit > 0
      d_sub <- d[keep, c("cellCode", "scientificName", "taxonKey", "year")]
      d_sub$occurrences <- hit[keep]
      d_sub
    })
    res_b <- suppressMessages(compute_impact_indicator(
      cube = bind_rows(drawn), impact_data = eicat_data,
      method = "mean_cum", ci_type = "none"
    ))
    vals <- res_b$impact
    mat[as.character(vals$year), b] <- vals$diversity_val
    if (b %% 10 == 0) cat(".")
  }
  cat(" done\n")

  # A year absent from a draw contributes a genuine zero, not a missing value.
  mat[is.na(mat)] <- 0

  tibble(
    year = as.integer(rownames(mat)),
    mean = rowMeans(mat),
    lo = apply(mat, 1, quantile, probs = 0.025, names = FALSE),
    hi = apply(mat, 1, quantile, probs = 0.975, names = FALSE)
  )
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

scale01 <- function(x, lo = min(x, na.rm = TRUE), hi = max(x, na.rm = TRUE)) {
  if (!is.finite(lo) || !is.finite(hi) || hi == lo) return(rep(NA_real_, length(x)))
  (x - lo) / (hi - lo)
}

# Fold-variation (max/min over non-zero values) as an effect-size summary.
fold_var <- function(x) {
  x <- x[!is.na(x) & x > 0]
  if (length(x) < 2) return(NA_real_)
  max(x) / min(x)
}

site_series <- list()
corr_tbl <- list()
coverage_log <- list()
fold_log <- list()

# --- Main loop -------------------------------------------------------------

for (skey in names(sites)) {
  s <- sites[[skey]]
  cat(sprintf("\n=== Site %d: %s ===\n", s$id, s$name))

  taxon_df <- read.delim(s$griis_taxon, stringsAsFactors = FALSE, sep = "\t", quote = "")
  profile_df <- read.delim(s$griis_profile, stringsAsFactors = FALSE, sep = "\t", quote = "")

  taxon_df <- taxon_df %>%
    mutate(clean_name = str_extract(scientificName, "^[^ ]+ [^ ]+"))

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

  # Identical to the size-based pipeline: a species counts as alien only if the
  # GRIIS species-profile join yields a non-missing isInvasive value.
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
  cat(sprintf("GRIIS alien: %d | matched to GIDIAS: %d\n",
              length(alien_species_list), length(eicat_species)))

  # Raw indicator on the full cube
  cat("Computing raw EICAT impact indicator...\n")
  res_imp <- compute_impact_indicator(
    cube = site_cube %>% rename(scientificName = species, cellCode = mgrscellcode,
                                taxonKey = specieskey),
    impact_data = eicat_data, method = "mean_cum", ci_type = "none"
  )
  impact_raw_df <- res_imp$impact %>% select(year, impact_raw = diversity_val)

  annual_effort <- site_cube %>%
    group_by(year) %>%
    summarise(total_occurrences = sum(occurrences), .groups = "drop")

  alien_rec <- site_cube %>%
    filter(species %in% eicat_species) %>%
    transmute(cellCode = mgrscellcode, scientificName = species,
              taxonKey = specieskey, year, occurrences)

  annual_alien <- alien_rec %>%
    group_by(year) %>%
    summarise(alien_records_eicat = sum(occurrences), .groups = "drop")

  # Record-level expansion, built once and shared by both rarefactions
  all_rare_years <- sort(annual_alien$year)
  year_rows <- lapply(all_rare_years, function(y) alien_rec[alien_rec$year == y, ])
  names(year_rows) <- as.character(all_rare_years)
  year_idx <- lapply(year_rows, function(d) rep.int(seq_len(nrow(d)), d$occurrences))

  # --- Coverage profile ---------------------------------------------------
  species_ab <- alien_rec %>%
    group_by(year, scientificName) %>%
    summarise(ab = sum(occurrences), .groups = "drop")

  cov_profile <- species_ab %>%
    group_by(year) %>%
    summarise(
      n_records = sum(ab),
      n_species = n(),
      f1 = sum(ab == 1),
      f2 = sum(ab == 2),
      coverage_observed = coverage_hat(ab),
      .groups = "drop"
    ) %>%
    arrange(year)

  cat("\n--- Per-year coverage profile (EICAT-matched alien records) ---\n")
  print(as.data.frame(cov_profile %>%
    mutate(coverage_observed = round(coverage_observed, 4))), row.names = FALSE)

  cov_ok <- cov_profile %>%
    filter(!is.na(coverage_observed), coverage_observed >= COVERAGE_FLOOR)
  cov_bad <- cov_profile %>%
    filter(is.na(coverage_observed) | coverage_observed < COVERAGE_FLOOR)
  zero_years <- setdiff(annual_effort$year, annual_alien$year)

  c_star <- if (nrow(cov_ok) > 0) min(cov_ok$coverage_observed) else NA_real_

  cat(sprintf("\nTarget coverage C* = %.4f | years retained = %d | excluded = %d\n",
              c_star, nrow(cov_ok), nrow(cov_bad) + length(zero_years)))
  if (nrow(cov_bad) > 0) {
    cat("  Excluded (coverage below floor / undefined): ",
        paste(sprintf("%d (C=%s)", cov_bad$year,
                      ifelse(is.na(cov_bad$coverage_observed), "NA",
                             sprintf("%.3f", cov_bad$coverage_observed))),
              collapse = ", "), "\n", sep = "")
  }
  if (length(zero_years) > 0) {
    cat("  Excluded (no EICAT-matched alien records): ",
        paste(sort(zero_years), collapse = ", "), "\n", sep = "")
  }

  # Solve for m in each retained year
  ab_by_year <- split(species_ab$ab, as.character(species_ab$year))
  cov_ok <- cov_ok %>%
    mutate(m_records_drawn = vapply(
      seq_len(n()),
      function(i) solve_m(ab_by_year[[as.character(year[i])]], c_star,
                          coverage_observed[i]),
      integer(1)
    ))

  cat(sprintf("m per year: min = %d, median = %.0f, max = %d (of n = %d..%d records)\n",
              min(cov_ok$m_records_drawn), median(cov_ok$m_records_drawn),
              max(cov_ok$m_records_drawn), min(cov_ok$n_records), max(cov_ok$n_records)))

  # --- Variant A(i): size-based rarefaction (reproduction check) ----------
  size_retained <- annual_alien %>% filter(alien_records_eicat >= N_MIN_FLOOR)
  n_size <- if (nrow(size_retained) > 0) min(size_retained$alien_records_eicat) else NA_integer_
  size_summary <- tibble(year = integer(0), mean = numeric(0),
                         lo = numeric(0), hi = numeric(0))
  if (!is.na(n_size) && nrow(size_retained) >= 3) {
    size_years <- as.character(sort(size_retained$year))
    m_size <- setNames(rep(as.integer(n_size), length(size_years)), size_years)
    size_summary <- run_rarefaction(year_rows, year_idx, m_size, eicat_data,
                                    seed = 20260729L + s$id, tag = "size-based")
  }

  # --- Variant A(ii): coverage-based rarefaction --------------------------
  cov_summary <- tibble(year = integer(0), mean = numeric(0),
                        lo = numeric(0), hi = numeric(0))
  if (nrow(cov_ok) >= 3) {
    cov_years <- as.character(cov_ok$year)
    m_cov <- setNames(as.integer(cov_ok$m_records_drawn), cov_years)
    cov_summary <- run_rarefaction(year_rows, year_idx, m_cov, eicat_data,
                                   seed = 20260731L + s$id, tag = "coverage-based")
  } else {
    cat("  Fewer than 3 years clear the coverage floor - coverage variant skipped.\n")
  }

  # --- Assemble -----------------------------------------------------------
  excl_reason <- bind_rows(
    cov_bad %>% transmute(
      year = as.integer(year),
      # as.character() keeps the column typed even when cov_bad has no rows
      year_excluded_reason = as.character(ifelse(
        is.na(coverage_observed),
        sprintf("excluded: coverage undefined (%d record(s))", n_records),
        sprintf("excluded: coverage %.3f below floor C = %.2f", coverage_observed, COVERAGE_FLOOR)))
    ),
    tibble(year = as.integer(zero_years),
           year_excluded_reason = rep("excluded: no EICAT-matched alien records",
                                      length(zero_years)))
  )

  series <- annual_effort %>%
    left_join(annual_alien, by = "year") %>%
    mutate(alien_records_eicat = coalesce(alien_records_eicat, 0L)) %>%
    left_join(impact_raw_df, by = "year") %>%
    left_join(cov_profile %>% select(year, coverage_observed), by = "year") %>%
    left_join(cov_ok %>% select(year, m_records_drawn), by = "year") %>%
    left_join(cov_summary %>% rename(impact_coverage_mean = mean,
                                     impact_coverage_lo = lo,
                                     impact_coverage_hi = hi), by = "year") %>%
    left_join(size_summary %>% rename(impact_size_mean = mean,
                                      impact_size_lo = lo,
                                      impact_size_hi = hi), by = "year") %>%
    left_join(excl_reason, by = "year") %>%
    mutate(
      target_coverage = ifelse(is.na(m_records_drawn), NA_real_, c_star),
      site_id = s$id, site_name = s$name
    ) %>%
    arrange(year)

  write_csv(
    series %>% select(year, total_occurrences, alien_records_eicat,
                      coverage_observed, target_coverage, m_records_drawn,
                      impact_coverage_mean, impact_coverage_lo, impact_coverage_hi,
                      year_excluded_reason, site_id, site_name),
    file.path(out_dir, sprintf("site_%d_invasive_impact_coverage.csv", s$id))
  )

  corr_tbl[[skey]] <- bind_rows(
    corr_row(s$id, s$name, "raw", series$impact_raw, series$total_occurrences),
    corr_row(s$id, s$name, "rarefied_size", series$impact_size_mean, series$total_occurrences),
    corr_row(s$id, s$name, "rarefied_coverage", series$impact_coverage_mean, series$total_occurrences)
  )

  # --- Sanity checks ------------------------------------------------------
  key <- as.character(s$id)
  got_raw <- corr_tbl[[skey]]$pearson_r_with_effort[1]
  got_size <- corr_tbl[[skey]]$pearson_r_with_effort[2]
  cat(sprintf("Sanity: raw r = %.4f (expect %.3f) | size-rarefied r = %.4f (expect %.3f)\n",
              got_raw, EXPECTED_RAW[[key]], got_size, EXPECTED_SIZE[[key]]))
  if (!is.finite(got_raw) || abs(got_raw - EXPECTED_RAW[[key]]) > 0.01) {
    stop(sprintf("SANITY CHECK FAILED (site %d): raw r = %.4f, expected %.3f.",
                 s$id, got_raw, EXPECTED_RAW[[key]]))
  }
  if (!is.finite(got_size) || abs(got_size - EXPECTED_SIZE[[key]]) > 0.01) {
    stop(sprintf("SANITY CHECK FAILED (site %d): size-rarefied r = %.4f, expected %.3f.",
                 s$id, got_size, EXPECTED_SIZE[[key]]))
  }

  coverage_log[[skey]] <- tibble(
    site_id = s$id, site_name = s$name,
    target_coverage = c_star,
    coverage_min_all_years = min(cov_profile$coverage_observed, na.rm = TRUE),
    coverage_max_all_years = max(cov_profile$coverage_observed, na.rm = TRUE),
    years_retained = nrow(cov_ok),
    years_excluded = nrow(cov_bad) + length(zero_years),
    m_min = min(cov_ok$m_records_drawn), m_median = median(cov_ok$m_records_drawn),
    m_max = max(cov_ok$m_records_drawn),
    size_rarefaction_n = as.integer(n_size),
    excluded_years = paste(sort(c(cov_bad$year, zero_years)), collapse = "; ")
  )

  fold_log[[skey]] <- tibble(
    site_id = s$id, site_name = s$name,
    fold_raw = fold_var(series$impact_raw),
    fold_size = fold_var(series$impact_size_mean),
    fold_coverage = fold_var(series$impact_coverage_mean)
  )

  site_series[[skey]] <- series
  rm(site_cube, alien_rec, year_rows, year_idx)
  gc(verbose = FALSE)
}

all_corr <- bind_rows(corr_tbl)
write_csv(
  all_corr %>% select(site_id, site_name, variant, pearson_r_with_effort,
                      spearman_rho, p_value, n_years),
  file.path(out_dir, "effort_standardisation_correlations_v2.csv")
)

all_series <- bind_rows(site_series)
cov_log <- bind_rows(coverage_log)
fold_tbl <- bind_rows(fold_log)

# Diagnostic: under coverage standardisation the draw size m varies by year, so
# check how much of the standardised series is simply tracking m.
m_diag <- all_series %>%
  filter(!is.na(m_records_drawn), !is.na(impact_coverage_mean)) %>%
  group_by(site_id, site_name) %>%
  summarise(
    r_m_vs_effort = cor(m_records_drawn, total_occurrences),
    r_impact_vs_m = cor(impact_coverage_mean, m_records_drawn),
    .groups = "drop"
  )

# --- Figure 10 (replaces the size-based version) ---------------------------

site_levels <- vapply(sites, function(x) x$name, character(1), USE.NAMES = FALSE)

series_long <- all_series %>%
  group_by(site_name) %>%
  mutate(
    `Raw impact (EICAT mean_cum)` = scale01(impact_raw),
    `Size-rarefied (constant n)` = scale01(impact_size_mean),
    `Coverage-rarefied (constant C*)` = scale01(impact_coverage_mean),
    `Total occurrences (effort)` = scale01(total_occurrences)
  ) %>%
  ungroup() %>%
  select(site_name, year, `Raw impact (EICAT mean_cum)`, `Size-rarefied (constant n)`,
         `Coverage-rarefied (constant C*)`, `Total occurrences (effort)`) %>%
  pivot_longer(-c(site_name, year), names_to = "series", values_to = "value_scaled") %>%
  # NA rows are kept deliberately: excluded years must show as gaps, not be bridged.
  mutate(site_name = factor(site_name, levels = site_levels))

ribbon_df <- all_series %>%
  group_by(site_name) %>%
  mutate(
    lo_scaled = scale01(impact_coverage_lo, min(impact_coverage_mean, na.rm = TRUE),
                        max(impact_coverage_mean, na.rm = TRUE)),
    hi_scaled = scale01(impact_coverage_hi, min(impact_coverage_mean, na.rm = TRUE),
                        max(impact_coverage_mean, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  filter(!is.na(impact_coverage_mean)) %>%
  select(site_name, year, lo_scaled, hi_scaled) %>%
  mutate(site_name = factor(site_name, levels = site_levels))

ann_df <- all_corr %>%
  select(site_name, variant, pearson_r_with_effort) %>%
  pivot_wider(names_from = variant, values_from = pearson_r_with_effort) %>%
  mutate(
    label = sprintf("raw r = %.2f  |  size-rarefied r = %.2f  |  coverage-rarefied r = %.2f",
                    raw, rarefied_size, rarefied_coverage),
    site_name = factor(site_name, levels = site_levels)
  )

series_levels <- c("Raw impact (EICAT mean_cum)", "Size-rarefied (constant n)",
                   "Coverage-rarefied (constant C*)", "Total occurrences (effort)")
series_colors <- setNames(c(ACCENT_RED, ACCENT_GREEN, ACCENT_BLUE, NEUTRAL_GREY), series_levels)
series_ltys <- setNames(c("solid", "solid", "solid", "dashed"), series_levels)
series_shapes <- setNames(c(16, 15, 17, NA), series_levels)
series_long$series <- factor(series_long$series, levels = series_levels)

# Title written from the observed result. Three things are tested: whether the
# standardised series still track effort, whether holding record count constant
# flattens the indicator, and whether the coverage series is just tracking its
# own varying draw size m.
r_of <- function(v) all_corr$pearson_r_with_effort[all_corr$variant == v]
n_sites <- length(unique(all_corr$site_id))
cov_still <- sum(r_of("rarefied_coverage") > 0.5, na.rm = TRUE)
cov_agrees <- sum(sign(r_of("rarefied_coverage")) == sign(r_of("rarefied_size")), na.rm = TRUE)
size_flat <- all(fold_tbl$fold_size < 0.1 * fold_tbl$fold_raw, na.rm = TRUE)
cov_m_driven <- all(abs(m_diag$r_impact_vs_m) > 0.95, na.rm = TRUE)

fig_title <- if (cov_still == 0 && (size_flat || cov_m_driven)) {
  paste0("Standardising for effort removes the sampling-volume signal - and leaves the ",
         "impact indicator with little independent variation")
} else if (cov_still == 0 && cov_agrees == n_sites) {
  "Coverage-based standardisation confirms it: invasive impact decouples from sampling volume"
} else if (cov_still == 0) {
  sprintf(paste0("Coverage-based standardisation decouples invasive impact from sampling volume, ",
                 "but disagrees with the size-based result at %d of %d sites"),
          n_sites - cov_agrees, n_sites)
} else {
  sprintf("Coverage-standardised invasive impact still tracks sampling volume at %d of %d sites",
          cov_still, n_sites)
}

fig_subtitle <- sprintf(
  paste0("Each series scaled 0-1 within its own site and metric so shapes are comparable. ",
         "Coverage-based rarefaction standardises to a common estimated sample coverage C* ",
         "(%s), drawing a year-specific number of records; size-based draws a constant n ",
         "(%s). B = %d draws; shaded band is the 95%% coverage-rarefaction interval."),
  paste(sprintf("%s: %.3f", cov_log$site_id, cov_log$target_coverage), collapse = ", "),
  paste(sprintf("%s: %d", cov_log$site_id, cov_log$size_rarefaction_n), collapse = ", "),
  B_DRAWS
)

p10 <- ggplot() +
  geom_ribbon(data = ribbon_df, aes(x = year, ymin = lo_scaled, ymax = hi_scaled),
              fill = ACCENT_BLUE, alpha = 0.15) +
  geom_line(data = series_long,
            aes(x = year, y = value_scaled, color = series, linetype = series),
            linewidth = 0.8) +
  geom_point(data = series_long %>% filter(series != "Total occurrences (effort)"),
             aes(x = year, y = value_scaled, color = series, shape = series), size = 1.3) +
  geom_text(data = ann_df, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.03, vjust = 1.4, size = 4.2, color = TEXT_COLOR,
            family = "sans", fontface = "italic") +
  facet_wrap(~site_name, ncol = 1) +
  scale_color_manual(values = series_colors, breaks = series_levels) +
  scale_linetype_manual(values = series_ltys, breaks = series_levels) +
  scale_shape_manual(values = series_shapes, breaks = series_levels, guide = "none") +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  # Clip rather than drop, so the rarefaction band stays continuous. The range is
  # matched to the coverage band (about -0.02 to 1.12) plus headroom for the panel
  # annotation; the old wider limits were sized for the size-based band and left
  # roughly 40% of each panel empty.
  coord_cartesian(ylim = c(-0.06, 1.38)) +
  labs(
    # Wrap widths are tied to the font sizes above - raising the sizes without
    # narrowing these would run the text off the canvas.
    title = str_wrap(fig_title, 70),
    subtitle = str_wrap(fig_subtitle, 88),
    x = "Year", y = "Scaled value (0-1 within series)",
    caption = str_wrap(paste0(
      "Source: GBIF occurrence cubes for three Ramsar sites; GRIIS national alien checklists; ",
      "GIDIAS EICAT assessments. Indicator: impIndicator::compute_impact_indicator(method = \"mean_cum\"). ",
      "Coverage estimator and interpolation follow Chao & Jost (2012). Annotations give the Pearson ",
      "correlation of each series with annual total occurrences."
    ), 125)
  ) +
  base_theme +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))

# Taller canvas at the same width: the header/legend/caption block grew with the
# larger type, so this gives the panels their height back. Embedding scales by
# width, so the on-page text size is unaffected.
ggsave(file.path(out_dir, "fig10_effort_standardised_impact.png"),
       plot = p10, width = 9.5, height = 12.0, dpi = 300, bg = "white")
ggsave(file.path(out_dir, "fig10_effort_standardised_impact.svg"),
       plot = p10, width = 9.5, height = 12.0, bg = "white")

fig_csv <- series_long %>%
  mutate(series_label = as.character(series), .keep = "unused") %>%
  left_join(ribbon_df %>% mutate(series_label = "Coverage-rarefied (constant C*)"),
            by = c("site_name", "year", "series_label")) %>%
  arrange(site_name, series_label, year)
write_csv(fig_csv, file.path(out_dir, "fig10_effort_standardised_impact.csv"))

# --- Console summary -------------------------------------------------------

cat("\n\n=============================================================\n")
cat("COVERAGE STANDARDISATION - CORRELATION WITH TOTAL OCCURRENCES\n")
cat("=============================================================\n")
print(as.data.frame(all_corr %>%
  mutate(pearson_r_with_effort = round(pearson_r_with_effort, 3),
         spearman_rho = round(spearman_rho, 3),
         p_value = signif(p_value, 3))), row.names = FALSE)

cat("\n--- Coverage settings actually used ---\n")
print(as.data.frame(cov_log %>%
  mutate(across(c(target_coverage, coverage_min_all_years, coverage_max_all_years),
                ~round(.x, 4))) %>%
  select(-excluded_years)), row.names = FALSE)

cat("\n--- Excluded years (coverage variant) ---\n")
for (i in seq_len(nrow(cov_log))) {
  cat(sprintf("  %d: %s\n", cov_log$site_id[i],
              ifelse(nzchar(cov_log$excluded_years[i]), cov_log$excluded_years[i], "none")))
}

cat("\n--- Diagnostic: is the coverage series just tracking its draw size m? ---\n")
print(as.data.frame(m_diag %>%
  mutate(across(starts_with("r_"), ~round(.x, 3)))), row.names = FALSE)

cat("\n--- Effect size: fold-variation (max/min) of each series ---\n")
print(as.data.frame(fold_tbl %>%
  mutate(across(starts_with("fold_"), ~round(.x, 1)))), row.names = FALSE)

cat("\n--- Files written ---\n")
cat("  report_figures/site_<id>_invasive_impact_coverage.csv (3 files)\n")
cat("  report_figures/effort_standardisation_correlations_v2.csv\n")
cat("  report_figures/fig10_effort_standardised_impact.png / .svg / .csv (replaced)\n")
cat("\nFigure title selected from the observed result:\n  ", fig_title, "\n", sep = "")
