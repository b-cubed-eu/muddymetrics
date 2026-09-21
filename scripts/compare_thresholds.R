#!/usr/bin/env Rscript
# Create comparison summary across thresholds

threshold_names <- c("liberal_0.10", "standard_0.25", "conservative_0.50")

cat("\n============================================================\n")
cat("       DATA SUFFICIENCY ANALYSIS - TROIA 2016 THRESHOLDS\n")
cat("                    (1990+ records)\n")
cat("============================================================\n\n")

# By continent comparison
cat("PASSING SITES BY CONTINENT\n")
cat("------------------------------------------------------------\n")
cat(sprintf("%-18s %10s %10s %10s\n", "Continent", "0.10/km²", "0.25/km²", "0.50/km²"))
cat("------------------------------------------------------------\n")

conts <- c("africa", "asia", "europe", "northamerica", "oceania", "southamerica")
for(cont in conts){
  vals <- c()
  for(tn in threshold_names){
    df <- read.csv(paste0("output/data_sufficiency/", tn, "/sites_passing.csv"))
    n <- sum(df$continent == cont)
    vals <- c(vals, n)
  }
  cat(sprintf("%-18s %10d %10d %10d\n", cont, vals[1], vals[2], vals[3]))
}

# Total
cat("------------------------------------------------------------\n")
vals <- c()
for(tn in threshold_names){
  df <- read.csv(paste0("output/data_sufficiency/", tn, "/sites_passing.csv"))
  vals <- c(vals, nrow(df))
}
cat(sprintf("%-18s %10d %10d %10d\n", "TOTAL", vals[1], vals[2], vals[3]))

cat("\n\n")

# Pass rate by continent
cat("PASS RATE (%) BY CONTINENT\n")
cat("------------------------------------------------------------\n")
cat(sprintf("%-18s %10s %10s %10s\n", "Continent", "0.10/km²", "0.25/km²", "0.50/km²"))
cat("------------------------------------------------------------\n")

for(cont in conts){
  vals <- c()
  for(tn in threshold_names){
    df <- read.csv(paste0("output/data_sufficiency/", tn, "/by_continent/summary.csv"))
    pct <- df$pass_rate[df$continent == cont]
    vals <- c(vals, ifelse(length(pct)>0, pct, 0))
  }
  cat(sprintf("%-18s %9.1f%% %9.1f%% %9.1f%%\n", cont, vals[1], vals[2], vals[3]))
}

cat("\n\n")

# Class comparison (top 10 classes in Europe)
cat("TOP 10 CLASSES BY PASS RATE (EUROPE)\n")
cat("------------------------------------------------------------\n")
cat(sprintf("%-25s %10s %10s %10s\n", "Class", "0.10/km²", "0.25/km²", "0.50/km²"))
cat("------------------------------------------------------------\n")

df_025 <- read.csv("output/data_sufficiency/standard_0.25/by_class/summary.csv")
euro_classes <- df_025[df_025$continent == "europe", ]
euro_classes <- euro_classes[order(-euro_classes$pass_rate), ]
top_classes <- head(euro_classes$class, 10)

for(cls in top_classes){
  vals <- c()
  for(tn in threshold_names){
    df <- read.csv(paste0("output/data_sufficiency/", tn, "/by_class/summary.csv"))
    sub <- df[df$class == cls & df$continent == "europe", ]
    if(nrow(sub) > 0){
      vals <- c(vals, sub$pass_rate)
    } else {
      vals <- c(vals, 0)
    }
  }
  cat(sprintf("%-25s %9.1f%% %9.1f%% %9.1f%%\n", cls, vals[1], vals[2], vals[3]))
}

cat("\n\n")

# ChaO2 summary
cat("CHAO2 COMPLETENESS (>=70% threshold)\n")
cat("------------------------------------------------------------\n")
chao2 <- read.csv("output/data_sufficiency/chao2_sac/chao2_sac_analysis.csv")
cont_chao <- aggregate(passes_chao2 ~ continent, data=chao2, FUN=sum)
cont_total <- aggregate(site_id ~ continent, data=chao2, FUN=length)
cont_merged <- merge(cont_chao, cont_total, by='continent')
cont_merged$pass_rate <- round(cont_merged$passes_chao2 / cont_merged$site_id * 100, 1)
print(cont_merged)

cat("\n============================================================\n")
cat("FILES SAVED IN:\n")
cat("  output/data_sufficiency/liberal_0.10/\n")
cat("  output/data_sufficiency/standard_0.25/\n")
cat("  output/data_sufficiency/conservative_0.50/\n")
cat("  output/data_sufficiency/chao2_sac/\n")
cat("============================================================\n")
