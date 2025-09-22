inputdir <- "inst/extdata/ramsar_site_data_100m_asia"
maindir <- "output/ramsar_metric_results_100m/asia"
shapefiledir <- "inst/extdata/ramsar_sites_wkt"
continent <- substring(inputdir, 36, nchar(inputdir))

total_occurrences <- calc_ramsar_indicator(
  "total_occ_ts", inputdir, maindir, shapefiledir, continent
)

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

