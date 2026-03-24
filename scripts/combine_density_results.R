library(dplyr)

outputdir <- "output/data_sufficiency"

files <- list.files(outputdir, pattern = "_density\\.csv$", full.names = TRUE)
results <- do.call(rbind, lapply(files, read.csv))

write.csv(results, file.path(outputdir, "all_sites_density.csv"), row.names = FALSE)

passing <- results %>% 
  filter(passes_threshold == TRUE) %>% 
  arrange(desc(density_records_km2))

write.csv(passing, file.path(outputdir, "sites_passing_density_threshold.csv"), row.names = FALSE)

cat("Total sites:", nrow(results), "\n")
cat("Passing:", nrow(passing), "\n")
cat("Pass rate:", round(nrow(passing)/nrow(results)*100, 1), "%\n")
cat("Mean density:", round(mean(results$density_records_km2, na.rm=TRUE), 4), "\n")
cat("Median density:", round(median(results$density_records_km2, na.rm=TRUE), 4), "\n")

# Show top 20 passing sites
cat("\nTop 20 sites by density:\n")
print(head(passing, 20))
