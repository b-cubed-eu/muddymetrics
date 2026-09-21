# Final Statistics Calculation
library(dplyr)

df <- read.csv("output/data_sufficiency/subgroup_analysis/sites_passing_all_precision.csv")

# Geographic continent mapping for all 116 countries
get_geographic_continent <- function(country) {
  case_when(
    country %in% c("Albania", "Andorra", "Austria", "Belarus", "Belgium", "Bosnia and Herzegovina", 
                   "Bulgaria", "Croatia", "Czechia", "Denmark", "Estonia", "Finland", "France", 
                   "Georgia", "Germany", "Hungary", "Iceland", "Ireland", "Italy", "Latvia", 
                   "Liechtenstein", "Lithuania", "Malta", "Monaco", "Montenegro", 
                   "Netherlands (Kingdom of the)", "North Macedonia", "Norway", "Poland", 
                   "Portugal", "Romania", "Serbia", "Slovakia", "Spain", "Sweden", "Switzerland", 
                   "Ukraine", "United Kingdom of Great Britain and Northern Ireland") ~ "Europe",
                   
    country %in% c("Algeria", "Benin", "Botswana", "Burkina Faso", "Cabo Verde", "Congo", 
                   "Cote dIvoire", "Democratic Republic of the Congo", "Egypt", "Equatorial Guinea", 
                   "Eswatini", "Gabon", "Ghana", "Guinea", "Guinea-Bissau", "Kenya", "Madagascar", 
                   "Malawi", "Mali", "Mauritania", "Morocco", "Mozambique", "Namibia", "Niger", 
                   "Nigeria", "Rwanda", "Sao Tome and Principe", "Senegal", "Seychelles", 
                   "South Africa", "Sudan", "Tunisia", "Zimbabwe") ~ "Africa",
                   
    country %in% c("Armenia", "Bhutan", "Cambodia", "China", "Cyprus", "India", "Indonesia", 
                   "Iran (Islamic Republic of)", "Iraq", "Japan", "Jordan", "Kazakhstan", 
                   "Kuwait", "Malaysia", "Myanmar", "Nepal", "Oman", "Philippines", 
                   "Republic of Korea", "Thailand", "Turkey", "United Arab Emirates", 
                   "Uzbekistan", "Vietnam") ~ "Asia",
                   
    country %in% c("Canada", "Dominican Republic", "El Salvador", "Guatemala", "Honduras", 
                   "Jamaica", "Mexico", "United States of America", "Panama") ~ "North America",
                   
    country %in% c("Argentina", "Brazil", "Chile", "Colombia", "Ecuador", 
                   "Peru", "Uruguay") ~ "South America",
                   
    country %in% c("Australia", "Fiji", "New Zealand", "Samoa", "Vanuatu") ~ "Oceania",
    
    TRUE ~ "Unknown"
  )
}

df$continent_geo <- get_geographic_continent(df$country)

# Unique passing sites by corrected continent
passing_sites_by_cont <- df %>%
  group_by(continent_geo) %>%
  summarize(
    passing_sites = n_distinct(site_id),
    .groups = "drop"
  )

# Get the denominators (total bounded sites) by corrected continent
wkt_base <- "inst/extdata/ramsar_sites_wkt"
wkt_files <- list.files(wkt_base, pattern = "\\.wkt$", recursive = TRUE, full.names = TRUE)

all_sites <- data.frame(
  wkt_path = wkt_files,
  stringsAsFactors = FALSE
) %>%
  mutate(
    country = basename(dirname(wkt_path))
  )

all_sites$continent_geo <- get_geographic_continent(all_sites$country)

denominators <- all_sites %>%
  group_by(continent_geo) %>%
  summarize(
    total_sites = n(),
    .groups = "drop"
  )

summary_table <- denominators %>%
  left_join(passing_sites_by_cont, by = "continent_geo") %>%
  mutate(
    passing_sites = ifelse(is.na(passing_sites), 0, passing_sites),
    pass_pct = (passing_sites / total_sites) * 100
  )

cat("Final Global Statistics for Presentation:\n")
cat(sprintf("%-15s | %-12s | %-12s | %-8s\n", "Continent", "Total Bounded", "Passing Sites", "Pass %"))
cat(rep("-", 55), "\n", sep="")

for (i in 1:nrow(summary_table)) {
  cat(sprintf("%-15s | %-12d | %-12d | %-7.2f%%\n", 
              summary_table$continent_geo[i], 
              summary_table$total_sites[i], 
              summary_table$passing_sites[i], 
              summary_table$pass_pct[i]))
}

cat(rep("-", 55), "\n", sep="")
total_bounded <- sum(summary_table$total_sites)
total_passing <- sum(summary_table$passing_sites)
cat(sprintf("%-15s | %-12d | %-12d | %-7.2f%%\n", "GLOBAL", total_bounded, total_passing, (total_passing/total_bounded)*100))

