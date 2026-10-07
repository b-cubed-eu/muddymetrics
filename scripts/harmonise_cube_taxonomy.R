#!/usr/bin/env Rscript
# Puts all site cubes on one taxonomy: the new GBIF backbone (Catalogue of Life
# based) that the 2024-2026 update uses. Rows up to 2023 carry the old backbone,
# where some species have other names (Sylvia curruca for Curruca curruca), so
# counting species by name across years splits them in two.
#
# Run after scripts/patch_site_cubes_update.R, in two steps:
#
#   1. Build the name table (needs internet; resumable):
#        Rscript scripts/harmonise_cube_taxonomy.R --build --limit 200   # first a small test
#        Rscript scripts/harmonise_cube_taxonomy.R --build               # then all names
#      Only names in the screened taxon groups are matched (birds, mammals,
#      herptiles, insects, plants and fish; --groups all matches every name).
#      Each old (kingdom, species) gets the accepted species of the new backbone:
#      names that occur in the update rows are taken from there, all others are
#      matched with the GBIF species-match API (v2, checklist
#      7ddf754f-d193-4cc9-b351-99906754a03b). API results are cached in
#      _local_data/taxonomy/api_matches.csv, so an interrupted run continues
#      where it stopped. The full table is _local_data/taxonomy/name_table.csv.
#
#   2. Apply it:
#        Rscript scripts/harmonise_cube_taxonomy.R --apply --sites site_67 --out-dir _local_data/taxonomy_check
#        Rscript scripts/harmonise_cube_taxonomy.R --apply              # all cubes, in place
#      Old-backbone rows get the new names and keys of all ranks (kingdom to
#      species); rows that then share cell, year and species are merged
#      (occurrences summed). Names that could not be matched keep their old
#      name and key. In-place runs copy each original to the backup folder first.
#
# Options:
#   --build | --apply | --report
#   --limit N        --build: match at most N names via the API in this run
#   --connections N  --build: parallel API requests [4]
#   --groups         --build: comma-separated taxon subsets whose names are matched,
#                    or "all" [SUBSETS in lib/sufficiency_criteria.R]
#   --sites          --apply: comma-separated site ids, or a file with one per line
#   --out-dir        --apply: write harmonised cubes here instead of in place
#   --backup-dir     --apply: [_local_data/cube_backup_taxonomy_<date>]
#   --dir            folder for the name table and cache [_local_data/taxonomy]

suppressPackageStartupMessages({ library(data.table); library(jsonlite); library(curl) })
source("scripts/lib/cube_update.R")
source("scripts/lib/taxonomy_match.R")
source("scripts/lib/sufficiency_criteria.R")   # SUBSETS, in_subsets()

arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
has <- function(flag) flag %in% commandArgs(trailingOnly = TRUE)
DIR <- arg("--dir", "_local_data/taxonomy")
ROOT <- "inst/extdata"
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
TABLE <- file.path(DIR, "name_table.csv")
CACHE <- file.path(DIR, "api_matches.csv")
cubes <- list_site_cubes(ROOT)

if (has("--build")) {
  cat(sprintf("Reading taxa from %d site cubes...\n", nrow(cubes)))
  acc <- list()
  for (f in cubes$file_path) {
    x <- fread(f, select = TAXON_COLS, colClasses = "character", showProgress = FALSE)
    x <- x[species != ""]
    acc[[f]] <- x[, .N, by = TAXON_COLS]
  }
  tx <- rbindlist(acc)[, .(n_rows = sum(N)), by = TAXON_COLS]
  tx[, old := is_old_key(specieskey)]
  old <- tx[old == TRUE, .(n_rows = sum(n_rows)), by = .(kingdom, phylum, class, order, species)]
  old <- old[order(-n_rows)][!duplicated(old[order(-n_rows)], by = c("kingdom", "species"))]
  groups <- arg("--groups", paste(SUBSETS, collapse = ","))
  if (groups != "all") {
    g <- strsplit(groups, ",")[[1]]
    n0 <- nrow(old)
    old <- old[in_subsets(old, g)]
    cat(sprintf("Matching names in %s only: %d of %d old-backbone names.\n", groups, nrow(old), n0))
  }
  newr <- tx[old == FALSE][order(-n_rows)][!duplicated(tx[old == FALSE][order(-n_rows)], by = c("kingdom", "species"))]
  cat(sprintf("%d old-backbone names (%s rows), %d new-backbone names.\n", nrow(old),
              format(sum(old$n_rows), big.mark = ","), nrow(newr)))
  lim <- as.numeric(arg("--limit", "Inf"))
  nt <- build_name_table(old, newr, CACHE, limit = lim,
                         fetch = function(u) fetch_urls(u, connections = as.integer(arg("--connections", "4"))))
  nt <- merge(nt, old[, .(old_kingdom = kingdom, old_species = species, n_rows)], by = c("old_kingdom", "old_species"), all.x = TRUE)
  fwrite(nt, TABLE)
  left <- nrow(old) - nrow(nt)
  cat(sprintf("Name table: %s (%d names; %d still to match via the API).\n", TABLE, nrow(nt), max(0, left)))
  print(nt[, .(names = .N, rows = sum(as.numeric(n_rows), na.rm = TRUE)), by = status])
  ch <- nt[status %in% c("accepted", "synonym") & species != old_species]
  if (nrow(ch)) { cat("Examples of names that change:\n"); print(head(ch[order(-as.numeric(n_rows)), .(old_species, species, class, status, match_type)], 15)) }
  quit(save = "no")
}

if (has("--report")) {
  nt <- fread(TABLE, colClasses = "character")
  print(nt[, .(names = .N, rows = sum(as.numeric(n_rows), na.rm = TRUE)), by = status])
  print(nt[status %in% c("accepted", "synonym", "in update data"), .(renamed = sum(species != old_species), lumped_into = sum(duplicated(paste(kingdom, species))))])
  quit(save = "no")
}

if (has("--apply")) {
  if (!file.exists(TABLE)) stop("No name table yet: run with --build first.")
  nt <- fread(TABLE, colClasses = "character")
  SITES <- arg("--sites"); OUT_DIR <- arg("--out-dir")
  BACKUP <- arg("--backup-dir", file.path("_local_data", paste0("cube_backup_taxonomy_", format(Sys.Date(), "%Y%m%d"))))
  todo <- cubes
  if (!is.null(SITES)) {
    w <- if (file.exists(SITES)) trimws(readLines(SITES)) else trimws(strsplit(SITES, ",")[[1]])
    todo <- todo[todo$site_id %in% w, ]
  }
  cat(sprintf("Harmonising %d site cubes%s.\n", nrow(todo), if (is.null(OUT_DIR)) " in place" else paste(" into", OUT_DIR)))
  summ <- list()
  for (i in seq_len(nrow(todo))) {
    f <- todo$file_path[i]
    x <- fread(f, colClasses = "character", showProgress = FALSE)
    n_old <- sum(is_old_key(x$specieskey))
    if (!n_old) { summ[[i]] <- data.table(site_id = todo$site_id[i], rows = nrow(x), old_rows = 0L, renamed = 0L, merged = 0L, unmatched_rows = 0L); next }
    h <- harmonise_rows(x, nt)
    unmatched <- sum(is_old_key(h$rows$specieskey))
    dest <- if (is.null(OUT_DIR)) f else file.path(OUT_DIR, sub(paste0("^", ROOT, "/?"), "", f))
    dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
    if (is.null(OUT_DIR)) {
      b <- file.path(BACKUP, sub(paste0("^", ROOT, "/?"), "", f))
      dir.create(dirname(b), showWarnings = FALSE, recursive = TRUE)
      if (!file.exists(b) && !file.copy(f, b, copy.date = TRUE)) stop("Backup failed for ", f)
    }
    tmpf <- paste0(dest, ".tmp"); fwrite(h$rows, tmpf)
    if (file.exists(dest)) file.remove(dest)
    if (!file.rename(tmpf, dest)) stop("Could not write ", dest)
    summ[[i]] <- data.table(site_id = todo$site_id[i], rows = nrow(x), old_rows = n_old, renamed = h$n_renamed,
                            merged = h$n_merged, unmatched_rows = unmatched)
    if (i %% 50 == 0) cat(sprintf("  %d of %d\n", i, nrow(todo)))
  }
  summ <- rbindlist(summ)
  fwrite(summ, file.path(DIR, "apply_summary.csv"))
  cat(sprintf("Done: %s old-backbone rows, %s renamed, %s merged, %s left with old names (unmatched).\n",
              format(sum(summ$old_rows), big.mark = ","), format(sum(summ$renamed), big.mark = ","),
              format(sum(summ$merged), big.mark = ","), format(sum(summ$unmatched_rows), big.mark = ",")))
  cat("Summary:", file.path(DIR, "apply_summary.csv"), "\n")
  quit(save = "no")
}

cat("Choose --build, --apply or --report (see the header of this script).\n")
