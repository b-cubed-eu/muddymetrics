inputdir <- "inst/extdata/ramsar_site_data_100m_asia"
maindir <- "output/ramsar_metric_results_100m/asia"
shapefiledir <- "inst/extdata/ramsar_sites_wkt"
continent <- substring(inputdir, 36, nchar(inputdir))

total_occurrences <- calc_ramsar_indicator(
  "total_occ_ts", inputdir, maindir, shapefiledir, continent
)

saveRDS(total_occurrences, file = paste0(maindir, "/asia_total_occurrences.rds"))

total_occurrences_map <- calc_ramsar_indicator(
  "total_occ_map", inputdir, maindir, shapefiledir, continent
)

saveRDS(total_occurrences_map, file = paste0(maindir, "/asia_total_occurrences_map.rds"))

observed_richness <- calc_ramsar_indicator(
  "total_occ_ts", inputdir, maindir, shapefiledir, continent
)

saveRDS(observed_richness, file = paste0(maindir, "/asia_observed_richness.rds"))

observed_richness_map <- calc_ramsar_indicator(
  "obs_richness_map", inputdir, maindir, shapefiledir, continent
)

saveRDS(observed_richness_map, file = paste0(maindir, "/asia_observed_richness_map.rds"))

cumulative_richness <- calc_ramsar_indicator(
  "cum_richness_ts", inputdir, maindir, shapefiledir, continent
)

saveRDS(cumulative_richness, file = paste0(maindir, "/asia_cumulative_richness.rds"))

cumulative_richness_map <- calc_ramsar_indicator(
  "cum_richness_map", inputdir, maindir, shapefiledir, continent
)

saveRDS(cumulative_richness_map, file = paste0(maindir, "/asia_cumulative_richness_map.rds"))




occ_per_year <- total_occurrences$values
mean_occ_per_year <- unlist(total_occurrences$mean)

howmany <- length(mean_occ_per_year > 100)

num_years2 <- unlist(num_years)
mean_num_years <- mean(num_years2, na.rm = TRUE)
howmany_numyears <- length(num_years2[num_years2 > 10])

num_years3 <- num_years2[names(num_years2) %in% names(enough_data2)]
temp_howmany <- vector()
for (i in seq_along(num_years3)) temp_howmany[i] <- (num_years3[[i]] > 10) && (enough_data2[[i]] > 100)
more_than_10_years_and_more_than_100_occperyear <- temp_howmany
saveRDS(more_than_10_years_and_more_than_100_occperyear, file = paste0(maindir, "/asia_more_than_10_years_and_more_than_100_occpearyear.rds"))

