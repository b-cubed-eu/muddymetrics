#!/usr/bin/env Rscript
library(readr)
library(dplyr)

density <- read_csv('output/data_sufficiency/density_by_class_1990.csv')

result <- density %>%
  group_by(continent, class) %>%
  summarise(total=n(), passing=sum(passes, na.rm=T), pct=round(passing/total*100,1), .groups='drop')

write.csv(result, 'output/data_sufficiency/class_continent_pass_rates.csv', row.names=F)

# Also create a summary by continent
cont_summary <- density %>%
  group_by(continent) %>%
  summarise(
    total_sites = n_distinct(site_id),
    total_combos = n(),
    passing_combos = sum(passes, na.rm=T),
    pass_rate = round(passing_combos/total_combos*100,1)
  )

write.csv(cont_summary, 'output/data_sufficiency/continent_summary.csv', row.names=F)

cat("Class-continent pass rates saved to output/data_sufficiency/\n")
