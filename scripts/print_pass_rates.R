#!/usr/bin/env Rscript
library(readr)
library(dplyr)
library(tidyr)

density <- read_csv('output/data_sufficiency/density_by_class_1990.csv')

# Get major classes (those with >=50 site-class combinations)
class_totals <- density %>% group_by(class) %>% summarise(n=n()) %>% filter(n >= 50)

# Calculate pass rate by class and continent for major classes
result <- density %>%
  filter(class %in% class_totals$class) %>%
  group_by(continent, class) %>%
  summarise(
    total = n(),
    passing = sum(passes, na.rm=T),
    pct = round(passing/total*100, 1),
    .groups = 'drop'
  )

# Pivot to wide
wide <- result %>% 
  pivot_wider(id_cols = class, names_from = continent, values_from = pct, values_fill = 0)

# Reorder columns and sort
wide <- wide[, c("class", "africa", "asia", "europe", "northamerica", "oceania", "southamerica")]
wide <- wide[order(-wide$europe), ]

# Print nicely
cat("\nPass rate (%) by class and continent (classes with >=50 site-class combos)\n")
cat("================================================================================\n")
cat(sprintf("%-25s %8s %8s %8s %12s %8s %12s\n", "Class", "Africa", "Asia", "Europe", "N. America", "Oceania", "S. America"))
cat("--------------------------------------------------------------------------------\n")

for(i in 1:nrow(wide)){
  row <- wide[i,]
  cat(sprintf("%-25s %7.1f%% %7.1f%% %7.1f%% %11.1f%% %7.1f%% %11.1f%%\n",
    row$class, row$africa, row$asia, row$europe, row$northamerica, row$oceania, row$southamerica))
}

cat("================================================================================\n\n")

# Also show overall by continent
cat("Overall pass rate by continent:\n")
cont <- density %>% group_by(continent) %>% summarise(pct=round(sum(passes)/n()*100,1), n=n())
print(cont)
