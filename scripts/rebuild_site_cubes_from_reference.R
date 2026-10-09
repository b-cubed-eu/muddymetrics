#!/usr/bin/env Rscript
# Rebuilds the rows before 2024 of the site cubes that lack many of them, from
# the rows collected by the reference step of scripts/check_site_cubes.R.
#
# The cube check of October 2026 found site cubes whose older rows are
# incomplete: cubes split on 16-17 June 2026 (reading stopped early, or the
# rows were taken from the wrong continent's cube), cubes of sites whose
# outlines are invalid in planar geometry, and rows filed under a neighbouring
# continent. The reference step has already read every line of the sources and
# kept the rows of every site, so the older rows are taken from there.
#
# For each site, the reference rows before --from-year are matched with the
# cube's rows (on every column except the higher-taxon counts). A cube is
# rebuilt when at least --min-rows of these reference rows, and at least
# --min-share of them, are missing from it. A site with reference rows but no
# cube gets a new cube. The new cube holds the reference rows before
# --from-year (all 25 columns as in the continental cubes) followed by the
# cube's own rows from --from-year on (for a new cube, the reference rows
# from --from-year on).
#
# Cubes filed under a continent folder other than the one whose continental
# cube holds most of the site's older rows (for example the Corsican site
# Etang de Biguglia under ramsar_site_data_100m_southamerica/France) are
# written to that continent's folder, same country folder and file name; the
# install step then moves the old file to the backup folder. Sites without
# older rows, or with a tie between continents, stay where they are.
# --no-refile leaves every cube where it is.
#
# The rebuilt cubes are written to --out-dir (same relative paths) and put
# into place with scripts/patch_site_cubes_update.R --install-from, which
# backs up and verifies each original. Run the cube check again afterwards.
#
# Usage (from the repository root):
#   Rscript scripts/rebuild_site_cubes_from_reference.R --dry-run      # list what would be rebuilt
#   Rscript scripts/rebuild_site_cubes_from_reference.R --out-dir _local_data/cube_rebuild
#
# Options:
#   --check-dir   output folder of the cube check [_local_data/cube_check]
#   --root        data folder [inst/extdata]
#   --out-dir     where the rebuilt cubes are written [_local_data/cube_rebuild]
#   --work-dir    per-bucket summaries and rebuild_summary.csv [_local_data/cube_rebuild_work]
#   --from-year   first year kept from the cubes [2024]
#   --min-rows    rebuild when at least this many older reference rows are missing [10]
#   --min-share   ... and at least this share of them [0.001]
#   --sites       only these site ids (comma-separated or a file)
#   --cores       parallel processes [4]
#   --force       rebuild every cube whose older rows differ from the reference at all
#                 (missing or extra rows), whatever --min-rows and --min-share say;
#                 meant with --sites
#   --no-refile   do not move cubes to the continent folder their rows come from
#   --dry-run     only write the summary of what would be done

suppressPackageStartupMessages({ library(data.table); library(sf); library(mgrs); library(stringi); library(parallel) })
source("scripts/lib/sufficiency_criteria.R")
source("scripts/lib/cube_check.R")
source("scripts/lib/cube_rebuild.R")

arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
has <- function(flag) flag %in% commandArgs(trailingOnly = TRUE)
CHECK <- arg("--check-dir", "_local_data/cube_check")
ROOT <- sub("/+$", "", arg("--root", "inst/extdata"))
OUT_DIR <- arg("--out-dir", "_local_data/cube_rebuild")
WORK <- arg("--work-dir", "_local_data/cube_rebuild_work")
FROM <- as.integer(arg("--from-year", "2024"))
RB_MIN_ROWS <- as.numeric(arg("--min-rows", "10"))
RB_MIN_SHARE <- as.numeric(arg("--min-share", "0.001"))
SITES <- arg("--sites")
CORES <- as.integer(arg("--cores", "4"))
DRY <- has("--dry-run")
RB_REFILE <- !has("--no-refile")
RB_FORCE <- has("--force")
REF_DIR <- file.path(CHECK, "reference")
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

if (!file.exists(file.path(REF_DIR, "sources_parts.csv"))) stop("No finished reference step in ", REF_DIR)
rs <- fread(file.path(REF_DIR, "sources_parts.csv"))
if (any(rs$from_year != FROM)) stop("The reference was built with --from-year ", paste(unique(rs$from_year), collapse = ","))
src <- fread(file.path(REF_DIR, "sources.csv"))
if (!all(src$read_completely)) stop("The reference step did not read every source completely; see ", file.path(REF_DIR, "sources.csv"))
NB <- unique(rs$n_buckets); if (length(NB) != 1) stop("The reference parts use different bucket counts")
if (normalizePath(OUT_DIR, mustWork = FALSE) == normalizePath(ROOT, mustWork = FALSE)) stop("--out-dir cannot be the data folder")
WORK <- if (DRY) file.path(WORK, "dry_run") else WORK
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
# Finished buckets are reused when the script is run again, so only with the same settings.
settings <- data.table(setting = c("check_dir", "root", "out_dir", "from_year", "min_rows", "min_share", "sites", "refile", "force"),
                       value = c(CHECK, ROOT, OUT_DIR, FROM, RB_MIN_ROWS, RB_MIN_SHARE, if (is.null(SITES)) "all" else SITES, RB_REFILE, RB_FORCE))
sf_ <- file.path(WORK, "settings.csv")
if (file.exists(sf_)) {
  old <- fread(sf_, colClasses = "character")
  if (!identical(old$value, as.character(settings$value))) stop(WORK, " holds results made with other settings; use another --work-dir")
} else fwrite(settings, sf_)

outlines <- readRDS(file.path(CHECK, "outlines.rds"))
cubes <- chk_list_cubes(ROOT)
cubes[, rel := substring(file, nchar(ROOT) + 2)]
wkt <- list.files(file.path(ROOT, "ramsar_sites_wkt"), pattern = "\\.wkt$", recursive = TRUE, full.names = TRUE)
ids <- unique(c(outlines$site_id, cubes$site_id[!is.na(cubes$site_id)]))
if (!is.null(SITES)) ids <- intersect(ids, if (file.exists(SITES)) trimws(readLines(SITES)) else trimws(strsplit(SITES, ",")[[1]]))
buckets <- sort(unique(chk_bucket_of(ids, NB)))
say("%d sites in %d buckets; rebuild when >= %g older rows and >= %.2f%% of them are missing; %s",
    length(ids), length(buckets), RB_MIN_ROWS, 100 * RB_MIN_SHARE, if (DRY) "dry run" else paste("writing to", OUT_DIR))

run_bucket <- function(b) {
  done <- file.path(WORK, sprintf("bucket_%03d.csv", b))
  if (file.exists(done)) return(b)
  ref_all <- rb_read_reference_bucket(REF_DIR, b)
  res <- lapply(ids[chk_bucket_of(ids, NB) == b], function(s) {
    files <- cubes$file[cubes$site_id %in% s]
    ref <- ref_all[site_id == s]
    row <- data.table(site_id = s, n_cube_files = length(files), cube = if (length(files)) paste(cubes$rel[cubes$site_id %in% s], collapse = " | ") else "",
                      rows_cube = 0L, rows_ref = nrow(ref), rows_before_cube = 0L, rows_before_ref = 0L, missing_before = 0L,
                      extra_before = 0L, rows_from_kept = 0L, rows_new = 0L, dup_keys_before = 0L,
                      continent_by_rows = "", moved_from = "",
                      sources = paste(sprintf("%s:%d", names(table(ref$source)), as.integer(table(ref$source))), collapse = "|"),
                      action = "kept")
    if (length(files) > 1) { row[, action := "skipped: several cube files"]; return(row) }
    cube <- if (length(files)) chk_read_cube(files) else NULL
    mb <- rb_missing_before(if (is.null(cube)) ref[0, CHK_COLS, with = FALSE] else cube, ref, FROM)
    n_cube <- if (is.null(cube)) 0L else nrow(cube)   # computed outside row[...], where `cube` would be the column
    row[, `:=`(rows_cube = n_cube, rows_before_cube = mb$rows_cube, rows_before_ref = mb$rows_ref,
               missing_before = mb$missing, extra_before = mb$extra)]
    rebuild <- if (is.null(cube)) nrow(ref) > 0 else if (RB_FORCE) mb$missing + mb$extra > 0 else
      mb$missing >= RB_MIN_ROWS && mb$missing >= RB_MIN_SHARE * mb$rows_ref
    cont_rows <- rb_continent_of(ref$source)
    row[, continent_by_rows := if (is.na(cont_rows)) "" else cont_rows]
    cur <- cubes[site_id %in% s]
    move <- !is.null(cube) && RB_REFILE && !is.na(cont_rows) && cont_rows != cur$folder
    if (!rebuild && !move) return(row)
    if (is.null(cube)) {
      w <- wkt[grepl(paste0("^", s, "_"), basename(wkt))]
      if (length(w) != 1) { row[, action := "skipped: no single outline file to name the new cube"]; return(row) }
      country <- basename(dirname(w))
      cont <- cont_rows
      if (is.na(cont)) {
        cc <- setdiff(unique(cubes$folder[cubes$country_dir == country]), "unknown")
        cont <- if (length(cc) == 1) cc else NA_character_
      }
      if (is.na(cont)) { row[, action := "skipped: continent folder not clear"]; return(row) }
      rel <- file.path(paste0("ramsar_site_data_100m_", cont), country, paste0(sub("\\.wkt$", "", basename(w)), "_data.csv"))
      row[, `:=`(cube = rel, action = "new cube")]
    } else {
      rel <- if (move) file.path(paste0("ramsar_site_data_100m_", cont_rows), cur$country_dir, basename(cur$file)) else cur$rel
      row[, `:=`(cube = rel, moved_from = if (move) cur$rel else "",
                 action = paste(c(if (rebuild) "rebuilt", if (move) "moved"), collapse = " and "))]
    }
    new <- if (rebuild) rb_rebuilt_rows(cube, ref, FROM) else cube
    yn <- suppressWarnings(as.numeric(new$year))
    row[, `:=`(rows_from_kept = sum(!is.na(yn) & yn >= FROM), rows_new = nrow(new),
               dup_keys_before = sum(duplicated(new[!is.na(yn) & yn < FROM, .(mgrscellcode, year, specieskey)])))]
    if (!DRY) {
      dest <- file.path(OUT_DIR, rel)
      dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
      tmp <- paste0(dest, ".tmp")
      fwrite(new, tmp)
      if (file.exists(dest)) file.remove(dest)
      if (!file.rename(tmp, dest)) stop("Could not write ", dest)
    }
    row
  })
  out <- rbindlist(res)
  fwrite(out, paste0(done, ".tmp")); file.rename(paste0(done, ".tmp"), done)
  b
}

cl <- makePSOCKcluster(CORES)
res <- tryCatch({
  # libraries first, then the settings, so that a name defined in a library cannot overwrite a setting
  clusterEvalQ(cl, {
    suppressPackageStartupMessages({ library(data.table); library(stringi) })
    source("scripts/lib/sufficiency_criteria.R"); source("scripts/lib/cube_check.R"); source("scripts/lib/cube_rebuild.R")
    NULL
  })
  clusterExport(cl, c("ROOT", "OUT_DIR", "WORK", "FROM", "RB_MIN_ROWS", "RB_MIN_SHARE", "RB_REFILE", "RB_FORCE", "DRY", "REF_DIR", "NB", "ids", "cubes", "wkt", "run_bucket"))
  parLapplyLB(cl, buckets, function(b) tryCatch(run_bucket(b), error = function(e) paste("bucket", b, ":", conditionMessage(e))), chunk.size = 1)
}, finally = stopCluster(cl))
errs <- unlist(res[vapply(res, is.character, TRUE)])
if (length(errs)) { print(errs); stop("Some buckets failed (listed above); finished buckets are kept, run again to continue.") }

summ <- rbindlist(lapply(file.path(WORK, sprintf("bucket_%03d.csv", buckets)), fread,
                         colClasses = list(character = c("site_id", "cube", "sources", "continent_by_rows", "moved_from"))))
summ <- summ[order(as.integer(sub("site_", "", site_id)))]
fwrite(summ, file.path(WORK, "rebuild_summary.csv"))
print(summ[, .(sites = .N, rows_before_ref = sum(rows_before_ref), missing_before = sum(missing_before), extra_before = sum(extra_before)), by = action])
say("Summary: %s", file.path(WORK, "rebuild_summary.csv"))
if (!DRY) say("Rebuilt cubes in %s; put them into place with: Rscript scripts/patch_site_cubes_update.R --install-from %s --backup-dir <folder> --work-dir %s",
              OUT_DIR, OUT_DIR, WORK)
