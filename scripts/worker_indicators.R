#!/usr/bin/env Rscript
# Computes the b3gbi time-series indicators and the MAGV reliability series
# for one site x taxon combination that passed the data-sufficiency screen.
#
# Usage (from the repository root):
#   Rscript scripts/worker_indicators.R <continent> <country> <site_id> <subset> <cube_path> <out_dir>
#
# Differences from worker_full_analysis.R (the D6.1 version):
#   - requires b3gbi >= 1.0.0 (0.9.x miscounts richness where species keys
#     are alphanumeric)
#   - species keys are harmonised by scientific name before process_cube(),
#     so a species recorded under two GBIF keys counts once
#   - herps include reptiles (Squamata, Testudines, Crocodylia, ...), using
#     the same subset definitions as the sufficiency screen
#   - the cube path is passed by the runner (no file-name matching)
#   - an indicator that fails with confidence intervals and falls back to
#     none is recorded in status.csv instead of silently
#   - MAGV is computed here from the total_occ series

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  cat("Usage: worker_indicators.R <continent> <country> <site_id> <subset> <cube_path> <out_dir>\n")
  quit(status = 1)
}
cont <- args[1]; country <- args[2]; site_id <- args[3]; subset <- args[4]
cube_path <- args[5]; out_dir <- args[6]

suppressPackageStartupMessages({
  library(b3gbi)
  library(data.table)
})
if (packageVersion("b3gbi") < "1.0.0") stop("b3gbi >= 1.0.0 is required; found ", packageVersion("b3gbi"))
source("scripts/lib/sufficiency_criteria.R")
source("scripts/lib/magv.R")

N_BOOT <- 1000
CI_TYPE <- "norm"
CELL_SIZE <- "auto"

# Faster, equivalent taxonomic distinctness (same result as
# b3gbi:::compute_tax_distinct_formula for any set of classifications:
# sum over pairs of shared ranks = sum over names of choose(count, 2)).
# Classifications that were not found (NA) are dropped.
fast_compute_tax_distinct_formula <- function(x, y) {
  tax_hier_temp <- y[names(y) %in% x$scientificName]
  tax_hier_temp <- tax_hier_temp[vapply(tax_hier_temp, is.data.frame, logical(1))]
  n_spec <- length(tax_hier_temp)
  if (n_spec < 3) return(NA)
  L <- max(vapply(tax_hier_temp, nrow, integer(1)))
  counts <- table(unlist(lapply(tax_hier_temp, function(df) unique(df$name))))
  shared_sum <- sum(choose(counts[counts > 1], 2))
  num_pairs <- n_spec * (n_spec - 1) / 2
  (L * num_pairs - shared_sum) / (L * num_pairs)
}
assignInNamespace("compute_tax_distinct_formula", fast_compute_tax_distinct_formula, ns = "b3gbi")

# Cache GBIF classifications so bootstrap replicates do not query GBIF again.
classification_cache <- list()
cached_classification <- function(x, ...) {
  missing <- setdiff(x, names(classification_cache))
  if (length(missing)) classification_cache[missing] <<- taxize::classification(missing, ...)
  classification_cache[x]
}
assignInNamespace("my_classification", cached_classification, ns = "b3gbi")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

df <- as.data.frame(data.table::fread(cube_path, showProgress = FALSE,
  colClasses = list(character = c("kingdom", "class", "species", "specieskey", "mgrscellcode"))))
df <- subset_taxon(filter_years(df), subset)
n_split <- names_with_multiple_keys(df)
df <- harmonise_species_keys(df)
if (nrow(df) < MIN_ROWS) quit(status = 0)

cube <- suppressMessages(suppressWarnings(b3gbi::process_cube(df,
  cols_cellCode = "mgrscellcode", cols_speciesKey = "specieskey", cols_scientificName = "species")))

indicators <- as.list(b3gbi::available_indicators)
status <- list()
for (ind in indicators) {
  if (is.null(ind$ts_wrapper) || is.na(ind$ts_wrapper)) next
  out_file <- file.path(out_dir, paste0(subset, "_", ind$indicator_class, "_ts.rds"))
  t0 <- Sys.time(); ci_used <- CI_TYPE; msg <- NA_character_
  res <- tryCatch(
    do.call(ind$ts_wrapper, list(data = cube, ci_type = CI_TYPE, num_bootstrap = N_BOOT, cell_size = CELL_SIZE)),
    error = function(e) { msg <<- conditionMessage(e); NULL })
  if (is.null(res)) {
    ci_used <- "none"
    res <- tryCatch(do.call(ind$ts_wrapper, list(data = cube, ci_type = "none", cell_size = CELL_SIZE)),
                    error = function(e) { msg <<- paste(msg, "| without CI:", conditionMessage(e)); NULL })
    if (is.null(res)) ci_used <- "failed"
  }
  if (!is.null(res)) saveRDS(res, out_file)
  status[[length(status) + 1]] <- data.frame(site_id = site_id, subset = subset,
    indicator = ind$indicator_class, ci = ci_used, message = msg,
    seconds = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  if (ind$indicator_class == "total_occ" && !is.null(res)) {
    occ <- as.data.frame(res$data)
    s <- magv_series(occ$year, occ$diversity_val)
    fwrite(cbind(data.frame(site_id = site_id, subset = subset), s),
           file.path(out_dir, paste0(subset, "_magv_by_year.csv")))
    fwrite(cbind(data.frame(site_id = site_id, subset = subset, continent = cont, country = country),
                 magv_summary(s)), file.path(out_dir, paste0(subset, "_magv_summary.csv")))
  }
}
st <- do.call(rbind, status)
st$b3gbi_version <- as.character(packageVersion("b3gbi"))
st$n_names_multi_key <- n_split
fwrite(st, file.path(out_dir, paste0(subset, "_status.csv")))
quit(status = 0)
