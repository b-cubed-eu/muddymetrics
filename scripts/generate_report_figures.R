#!/usr/bin/env Rscript
# generate_report_figures.R
# Generates publication-quality figures 4, 5, and 6 as PNG and SVG in report_figures/

library(ggplot2)
library(dplyr)
library(patchwork)
library(svglite)
library(ggrepel)
library(tools)

# 0. Setup
out_dir <- "report_figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Shared Styling Constants
TEXT_COLOR <- "#202124"
TITLE_COLOR <- "#1A365D" # Bold dark-slate
SUBTITLE_COLOR <- "#5F6368"
GRID_COLOR <- "#F1F3F4"
ACCENT_GREEN <- "#3B8A57" # Muted green for pass/increase
NEUTRAL_GREY <- "#9AA0A6" # Neutral grey for no change
ACCENT_RED <- "#C1584B"   # Muted red for decrease
SLATE_BLUE <- "#4A76A8"   # Slate blue for funnel bars

base_theme <- theme_minimal(base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 16, color = TITLE_COLOR, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 11, color = SUBTITLE_COLOR, margin = margin(b = 12)),
    plot.title.position = "plot",
    axis.title = element_text(face = "bold", size = 10, color = TEXT_COLOR),
    axis.text = element_text(size = 9, color = TEXT_COLOR),
    panel.grid.major = element_line(color = GRID_COLOR),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

# ==========================================
# 1. Figure 4 — The data-sufficiency funnel
# ==========================================
cat("--- Figure 4 ---\n")
funnel_csv_path <- "output/data_sufficiency/subgroup_analysis/all_tests_results_corrected.csv"
if (!file.exists(funnel_csv_path)) {
  stop("Input for Fig 4 not found!")
}

df_funnel <- read.csv(funnel_csv_path, stringsAsFactors = FALSE)
total_combos <- nrow(df_funnel)

# Compute cumulative survival
c0 <- total_combos
c1 <- sum(df_funnel$pass_density)
c2 <- sum(df_funnel$pass_density & df_funnel$pass_chao2)
c3 <- sum(df_funnel$pass_density & df_funnel$pass_chao2 & df_funnel$pass_sac_slope)
c4 <- sum(df_funnel$pass_density & df_funnel$pass_chao2 & df_funnel$pass_sac_slope & df_funnel$pass_temporal)
c5 <- sum(df_funnel$pass_density & df_funnel$pass_chao2 & df_funnel$pass_sac_slope & df_funnel$pass_temporal & df_funnel$pass_precision)

fig4_data <- data.frame(
  StageIndex = 6:1,
  Stage = factor(c("Start (all combinations)", 
                   "Observation density (>= 0.25 occ/km2)", 
                   "Chao2 completeness (>= 0.70)", 
                   "SAC slope (<= 0.10)", 
                   "Temporal decoupling (Jaccard elbow)", 
                   "Spatial precision ratio (<= 0.10)"),
                 levels = rev(c("Start (all combinations)", 
                            "Observation density (>= 0.25 occ/km2)", 
                            "Chao2 completeness (>= 0.70)", 
                            "SAC slope (<= 0.10)", 
                            "Temporal decoupling (Jaccard elbow)", 
                            "Spatial precision ratio (<= 0.10)"))),
  Count = c(c0, c1, c2, c3, c4, c5),
  stringsAsFactors = FALSE
)
fig4_data$Pct <- fig4_data$Count / total_combos * 100

# Write CSV
write.csv(fig4_data, file.path(out_dir, "fig4_sufficiency_funnel.csv"), row.names = FALSE)
cat(sprintf("Saved fig4 CSV. Total combos: %d, Passing: %d\n", total_combos, c5))

# Plot Fig 4
p4 <- ggplot(fig4_data, aes(y = Stage, x = Count)) +
  geom_bar(stat = "identity", aes(fill = Stage == "Spatial precision ratio (<= 0.10)"), width = 0.6) +
  scale_fill_manual(values = c("FALSE" = SLATE_BLUE, "TRUE" = ACCENT_GREEN)) +
  geom_text(aes(label = sprintf("%s (%.2f%%)", format(Count, big.mark=","), Pct)), 
            hjust = -0.1, color = TEXT_COLOR, size = 3.5, fontface = "bold") +
  
  # Chao2 drop annotation (text left-aligned at tail)
  annotate("text", x = 4050, y = 4.25, label = "Chao2 completeness removes ~58% of combinations", 
           color = ACCENT_RED, fontface = "italic", size = 3.2, hjust = 0) +
  
  # SAC slope drop annotation (text left-aligned at tail)
  annotate("text", x = 1850, y = 3.25, label = "SAC-slope saturation removes ~82% of remainder", 
           color = ACCENT_RED, fontface = "italic", size = 3.2, hjust = 0) +
  
  # End annotation
  annotate("label", x = 4000, y = 1.3, 
           label = "77 passing combinations = 60 unique sites\n(3.7% of 1,638 mapped sites)", 
           fill = "#F8F9FA", color = TITLE_COLOR, size = 3.5, fontface = "bold", label.size = 0.2) +
  
  scale_x_continuous(limits = c(0, 9200), expand = c(0, 0), labels = scales::comma) +
  labs(
    title = "Most Ramsar sites fail open-data sufficiency on depth",
    subtitle = "Cumulative survival of 6,761 network-wide combinations across sequential filtering stages",
    x = "Surviving Combinations",
    y = ""
  ) +
  base_theme +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    plot.margin = margin(15, 30, 15, 15)
  )

ggsave(file.path(out_dir, "fig4_sufficiency_funnel.png"), plot = p4, width = 9.5, height = 5.5, dpi = 300)
ggsave(file.path(out_dir, "fig4_sufficiency_funnel.svg"), plot = p4, width = 9.5, height = 5.5)


# ==========================================
# 2. Figure 5 — The "success illusion"
# ==========================================
cat("--- Figure 5 ---\n")
trend_csv_path <- "output/part2_trend_direction.csv"
if (!file.exists(trend_csv_path)) {
  stop("Input for Fig 5 not found!")
}

df_trend <- read.csv(trend_csv_path, stringsAsFactors = FALSE)

# Explicit list of the 12 target indicators
target_indicators <- c(
  "total_occ", "occ_density", "cum_richness", "obs_richness",
  "hill0", "hill1", "hill2", "pielou_evenness", "williams_evenness",
  "ab_rarity", "occ_turnover", "tax_distinct"
)

# Aggregate globally
fig5_agg <- df_trend %>%
  filter(indicator %in% target_indicators) %>%
  group_by(indicator) %>%
  summarise(
    n_pos = sum(n_positive),
    n_neg = sum(n_negative),
    n_neu = sum(n_neutral),
    n_total = sum(n_combos),
    .groups = "drop"
  ) %>%
  mutate(
    pct_pos = n_pos / n_total * 100,
    pct_neg = n_neg / n_total * 100,
    pct_neu = n_neu / n_total * 100
  )

# Add group and DisplayName
fig5_agg <- fig5_agg %>%
  mutate(
    Group = factor(case_when(
      indicator %in% c("total_occ", "occ_density", "cum_richness", "obs_richness") ~ "Effort-Sensitive Indicators",
      TRUE ~ "Standardised & Structural Indicators"
    ), levels = c("Effort-Sensitive Indicators", "Standardised & Structural Indicators")),
    DisplayName = case_when(
      indicator == "total_occ" ~ "Total occurrences",
      indicator == "occ_density" ~ "Occurrence density",
      indicator == "cum_richness" ~ "Cumulative richness",
      indicator == "obs_richness" ~ "Observed richness",
      indicator == "hill0" ~ "Estimated richness (Hill q=0)",
      indicator == "hill1" ~ "Hill-Shannon (q=1)",
      indicator == "hill2" ~ "Hill-Simpson (q=2)",
      indicator == "pielou_evenness" ~ "Pielou evenness",
      indicator == "williams_evenness" ~ "Williams evenness",
      indicator == "ab_rarity" ~ "Abundance-based rarity",
      indicator == "occ_turnover" ~ "Occupancy turnover",
      indicator == "tax_distinct" ~ "Taxonomic distinctness",
      TRUE ~ indicator
    )
  )

# Set factor levels for correct vertical sorting
indicator_order <- c(
  "Total occurrences",
  "Occurrence density",
  "Cumulative richness",
  "Observed richness",
  "Estimated richness (Hill q=0)",
  "Hill-Shannon (q=1)",
  "Hill-Simpson (q=2)",
  "Pielou evenness",
  "Williams evenness",
  "Abundance-based rarity",
  "Occupancy turnover",
  "Taxonomic distinctness"
)
fig5_agg$DisplayName <- factor(fig5_agg$DisplayName, levels = rev(indicator_order))

# Write CSV
write.csv(fig5_agg, file.path(out_dir, "fig5_success_illusion.csv"), row.names = FALSE)

# Generate split neutral dataset for geom_col(position="stack")
fig5_plot_data <- bind_rows(
  fig5_agg %>% mutate(FillCategory = "Decreasing", Value = -pct_neg),
  fig5_agg %>% mutate(FillCategory = "Neutral_Neg", Value = -0.5 * pct_neu),
  fig5_agg %>% mutate(FillCategory = "Neutral_Pos", Value = 0.5 * pct_neu),
  fig5_agg %>% mutate(FillCategory = "Increasing", Value = pct_pos)
)

fig5_plot_data$FillCategory <- factor(
  fig5_plot_data$FillCategory,
  levels = c("Neutral_Neg", "Neutral_Pos", "Decreasing", "Increasing")
)

# Generate labels dataset
fig5_labels <- bind_rows(
  fig5_agg %>% mutate(x = -0.5 * pct_neu - 0.5 * pct_neg, label = ifelse(n_neg > 2, as.character(n_neg), "")),
  fig5_agg %>% mutate(x = 0, label = ifelse(n_neu > 0, as.character(n_neu), "")),
  fig5_agg %>% mutate(x = 0.5 * pct_neu + 0.5 * pct_pos, label = ifelse(n_pos > 0, as.character(n_pos), ""))
) %>%
  filter(label != "")

# Create Plot
p5 <- ggplot(fig5_plot_data, aes(x = Value, y = DisplayName, fill = FillCategory)) +
  geom_col(width = 0.55, position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(
    name = NULL,
    breaks = c("Decreasing", "Neutral_Pos", "Increasing"),
    labels = c("Decreasing (↓)", "No Change (—)", "Increasing (↑)"),
    values = c(
      "Decreasing" = ACCENT_RED,
      "Neutral_Neg" = NEUTRAL_GREY,
      "Neutral_Pos" = NEUTRAL_GREY,
      "Increasing" = ACCENT_GREEN
    )
  ) +
  
  # Add count labels
  geom_text(
    data = fig5_labels,
    aes(x = x, y = DisplayName, label = label),
    color = "white",
    size = 3.2,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  # Add a vertical center line
  geom_vline(xintercept = 0, linetype = "dashed", color = TEXT_COLOR, alpha = 0.5) +
  
  facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(
    limits = c(-100, 100),
    breaks = seq(-100, 100, 25),
    labels = function(x) paste0(abs(x), "%")
  ) +
  
  labs(
    title = "More looking, not more biodiversity: raw counts rise while standardised diversity stays flat",
    subtitle = "Global trend direction counts across the 77 data-sufficient combinations (linear slope p <= 0.1)",
    x = "Percentage of Combinations",
    y = ""
  ) +
  base_theme +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "bold", size = 10, color = TITLE_COLOR),
    strip.background = element_rect(fill = "#F1F3F4", color = NA),
    panel.spacing = unit(15, "pt"),
    plot.margin = margin(15, 100, 15, 15)
  )

ggsave(file.path(out_dir, "fig5_success_illusion.png"), plot = p5, width = 10.5, height = 8.0, dpi = 300)
ggsave(file.path(out_dir, "fig5_success_illusion.svg"), plot = p5, width = 10.5, height = 8.0)


# ==========================================
# 3. Figure 6 — Monitoring reliability (MAGV)
# ==========================================
cat("--- Figure 6 ---\n")
magv_csv_path <- "output/magv_summary.csv"
if (!file.exists(magv_csv_path)) {
  stop("Input for Fig 6 not found!")
}

df_magv <- read.csv(magv_csv_path, stringsAsFactors = FALSE)

# Generate CSV
write.csv(df_magv, file.path(out_dir, "fig6_monitoring_reliability.csv"), row.names = FALSE)

# Panel A: Lollipop plot of stable run lengths
df_panelA <- df_magv %>%
  arrange(desc(stable_run_length), desc(mean_reliability)) %>%
  mutate(
    Index = row_number(),
    Highlight = stable_run_length >= 5,
    Label = case_when(
      site_id == "site_395" & subset == "aves" ~ "UK site_395 (aves): 5 yr",
      site_id == "site_1252" & subset == "aves" ~ "NL site_1252 (aves): 6 yr",
      site_id == "site_1279" & subset == "aves" ~ "NL site_1279 (aves): 5 yr",
      TRUE ~ ""
    )
  )

p6a <- ggplot(df_panelA, aes(x = reorder(paste(site_id, subset, sep="-"), -stable_run_length), y = stable_run_length)) +
  geom_segment(aes(xend = reorder(paste(site_id, subset, sep="-"), -stable_run_length), yend = 0), color = "#E0E0E0", linewidth = 0.5) +
  geom_point(aes(color = Highlight, size = Highlight)) +
  scale_color_manual(values = c("FALSE" = NEUTRAL_GREY, "TRUE" = ACCENT_GREEN)) +
  scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 4)) +
  geom_hline(yintercept = 5, linetype = "dashed", color = ACCENT_RED, alpha = 0.8) +
  
  # Text labels for highlights using ggrepel
  geom_text_repel(
    aes(label = Label),
    color = TITLE_COLOR,
    size = 3.2,
    fontface = "bold",
    nudge_y = 1.2,
    box.padding = 0.6,
    point.padding = 0.4,
    segment.color = TITLE_COLOR,
    segment.size = 0.3,
    arrow = arrow(length = unit(0.01, "npc")),
    force = 2,
    max.overlaps = Inf
  ) +
  
  scale_y_continuous(limits = c(0, 8), breaks = 0:8, expand = c(0.1, 0.1)) +
  labs(
    title = "Panel A: Continuous stable monitoring window length",
    subtitle = "Length of continuous years with moving average variability (MAGV) reliability >= 0.80",
    x = "77 Data-Sufficient Combinations",
    y = "Stable Window Length (Years)"
  ) +
  base_theme +
  theme(
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )

# Panel B: Density & Rug of mean_reliability
p6b <- ggplot(df_magv %>% filter(!is.na(mean_reliability)), aes(x = mean_reliability)) +
  geom_density(fill = "#E8F0FE", color = SLATE_BLUE, linewidth = 1) +
  geom_rug(color = TITLE_COLOR, sides = "b", length = unit(0.08, "npc")) +
  
  # Median line
  geom_vline(xintercept = 0.4624, linetype = "dashed", color = TEXT_COLOR, alpha = 0.8) +
  annotate("text", x = 0.44, y = 1.8, label = "Median = 0.46", color = TEXT_COLOR, angle = 90, size = 3.5, fontface = "bold") +
  
  # Max value line
  geom_vline(xintercept = 0.7187, linetype = "dashed", color = SLATE_BLUE, alpha = 0.8) +
  annotate("text", x = 0.74, y = 1.8, label = "Max = 0.72", color = SLATE_BLUE, angle = 90, size = 3.5, fontface = "bold") +
  
  # Stability Threshold
  geom_vline(xintercept = 0.80, linetype = "solid", color = ACCENT_RED, size = 1) +
  annotate("text", x = 0.83, y = 2.0, label = "Stability Threshold (0.80)\n[No site reaches this]", 
           color = ACCENT_RED, size = 3.5, fontface = "bold", hjust = 0) +
  
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Panel B: Distribution of mean monitoring reliability",
    subtitle = "Mean MAGV reliability score across full time series for all 77 combinations",
    x = "Mean Monitoring Reliability",
    y = "Density"
  ) +
  base_theme +
  theme(
    panel.grid.major.y = element_blank()
  )

# Combine using patchwork
p6 <- p6a / p6b + 
  plot_annotation(
    title = "Even the data-sufficient sites rarely offer a stable window to read a trend",
    subtitle = "Only 3 of the 77 passing combinations sustain a stable window of >= 5 years, and none reaches 0.80 reliability",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, color = TITLE_COLOR, margin = margin(t=15, b = 4)),
      plot.subtitle = element_text(size = 11, color = SUBTITLE_COLOR, margin = margin(b = 12)),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

ggsave(file.path(out_dir, "fig6_monitoring_reliability.png"), plot = p6, width = 10.0, height = 8.5, dpi = 300)
ggsave(file.path(out_dir, "fig6_monitoring_reliability.svg"), plot = p6, width = 10.0, height = 8.5)

cat("\nDone! Figures 4, 5, and 6 successfully generated in report_figures/.\n")
