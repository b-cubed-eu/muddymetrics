#!/usr/bin/env Rscript
# Adds the 2024-2026 GBIF update to the site cubes it never reached.
#
# The June 2026 update (GBIF download 0074420-260519110011954) was written into
# 1,147 site cubes, but sites were skipped when their country has folders under
# more than one continent, when they were re-split for boundary changes on 16-17
# June, or when their name contains a dot. This script finds every site cube
# that still lacks the update and rebuilds its rows from 2024 onwards from the
# update file: years before 2024 are kept, rows from 2024 on are replaced.
#
# A site cube "lacks the update" when it has no rows from 2024 on, or when most
# of them carry old-backbone (numeric) species keys; rows from the update carry
# new-backbone keys (short codes such as "7CYVX"). Site lists come from the cube files and outline
# files themselves (no country-to-continent lookup), and site ids are read from
# file names without splitting on ".".
#
# Usage (from the repository root):
#   Rscript scripts/patch_site_cubes_update.R --dry-run          # list the sites only
#   Rscript scripts/patch_site_cubes_update.R --all              # rebuild 2024+ in every cube (recommended)
#   Rscript scripts/patch_site_cubes_update.R --all --update <main update>.csv,<fish download>.csv \
#     --update-groups all,fish --to-year 2025
#        # the same, adding the fish of a second download, years 2024-2025 only
#   Rscript scripts/patch_site_cubes_update.R                    # only the cubes the update never reached
#   Rscript scripts/patch_site_cubes_update.R --sites site_67,site_838 --out-dir _local_data/cube_patch_check
#        # check run: write patched copies elsewhere, leave the cubes untouched
#   Rscript scripts/patch_site_cubes_update.R --all --out-dir _local_data/cube_patch_all ...
#   Rscript scripts/patch_site_cubes_update.R --install-from _local_data/cube_patch_all --dry-run
#   Rscript scripts/patch_site_cubes_update.R --install-from _local_data/cube_patch_all
#        # patch every cube into a separate folder, check it, then put the cubes into place
#        # (each original is copied to the backup folder first; nothing is read from the updates)
#
# Options:
#   --update      update cube(s) (tab-separated, 25 columns); several files are given
#                 comma-separated, e.g. the main update followed by a download of
#                 groups it does not cover (fish)
#                 [inst/extdata/ramsar_site_data_100m_world_2024_2026/0074420-260519110011954.csv]
#   --update-groups  comma-separated, one per update file: the taxon subsets
#                 (lib/sufficiency_criteria.R) taken from that file, "all" for every row,
#                 several joined with "+" (e.g. aves+mammalia) [all for every file]
#   --from-year   first year taken from the update [2024]
#   --to-year     last year taken from the update; later rows, in the update or in
#                 the existing cube, are left out [no limit]
#   --all         rebuild every site cube, not only those that lack the update
#   --sites       comma-separated site ids, or a file with one per line: patch these
#                 sites whether or not they need it (for checks)
#   --out-dir     write patched cubes here (same relative paths) instead of in place
#   --backup-dir  in-place runs copy each original here first [_local_data/cube_backup_<date>]
#   --work-dir    temporary per-site files [_local_data/cube_patch_work]
#   --chunk-lines lines read per chunk [2000000]
#   --dry-run     list the sites that would be patched (or, with --install-from, the
#                 cubes that would be replaced) and stop
#   --install-from  put the cubes of an earlier --out-dir run into place (see above)
#
# Outputs (in the work dir): patch_targets.csv (sites and why), patch_summary.csv
# (rows before and after, per site, and how many came from the second and later
# update files), and for already-updated sites whether the rows rebuilt from the
# first update file are identical to the existing ones.
#
# The update also dropped all reptiles (class Reptilia in the new backbone) from
# the cubes it did reach, apparently through a class list with the old-backbone
# reptile classes (Squamata, Testudines, Crocodylia). Rebuilding every cube with
# --all restores them; for cubes that already had the update, the summary checks
# that all other rebuilt rows are identical to the existing ones.
#
# Sites whose cube sits in ramsar_site_data_100m_unknown (site 1979, misplaced by
# the update because of the dot in its name) are rebuilt into their country's
# continent folder; the misplaced file is moved to the backup folder.

suppressPackageStartupMessages({ library(data.table); library(sf); library(mgrs); library(stringi) })
source("scripts/lib/cube_update.R")
source("scripts/lib/sufficiency_criteria.R")   # ALL_SUBSETS, in_subsets(), for --update-groups

arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
has <- function(flag) flag %in% commandArgs(trailingOnly = TRUE)
UPDATE <- trimws(strsplit(arg("--update", "inst/extdata/ramsar_site_data_100m_world_2024_2026/0074420-260519110011954.csv"), ",")[[1]])
FROM <- as.integer(arg("--from-year", "2024"))
TO <- as.numeric(arg("--to-year", "Inf"))
GROUPS <- trimws(strsplit(arg("--update-groups", paste(rep("all", length(UPDATE)), collapse = ",")), ",")[[1]])
if (length(GROUPS) != length(UPDATE)) stop("--update-groups needs one entry per update file")
GROUPS <- lapply(strsplit(GROUPS, "+", fixed = TRUE), trimws)
bad <- setdiff(unlist(GROUPS), c("all", ALL_SUBSETS))
if (length(bad)) stop("Unknown taxon subset in --update-groups: ", paste(bad, collapse = ", "))
SITES <- arg("--sites")
OUT_DIR <- arg("--out-dir")
BACKUP <- arg("--backup-dir", file.path("_local_data", paste0("cube_backup_", format(Sys.Date(), "%Y%m%d"))))
WORK <- arg("--work-dir", "_local_data/cube_patch_work")
CHUNK <- as.numeric(arg("--chunk-lines", "2e6"))
ROOT <- "inst/extdata"; WKT_ROOT <- file.path(ROOT, "ramsar_sites_wkt")
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)

# --- Installing the cubes of an earlier --out-dir run --------------------------------
INSTALL <- arg("--install-from")
if (!is.null(INSTALL)) {
  log <- install_patched(INSTALL, ROOT, BACKUP, dry_run = has("--dry-run"))
  # time-stamped, so that a later install with the same work dir does not overwrite the log
  f <- file.path(WORK, sprintf("%s_%s.csv", if (has("--dry-run")) "install_plan" else "install_log", format(Sys.time(), "%Y%m%d_%H%M%S")))
  fwrite(log, f)
  print(log[, .N, by = action])
  if (!has("--dry-run")) cat(sprintf("Cubes in place: %d, all verified: %s. Originals in %s. Log: %s\n",
                                     nrow(log), all(log$verified), BACKUP, f))
  if (any(nzchar(log$strays))) { cat("Other files of the same sites, moved to the backup folder:\n"); print(log[nzchar(strays), .(site_id, strays)]) }
  quit(save = "no")
}

if (any(!file.exists(UPDATE))) stop("Update file not found: ", paste(UPDATE[!file.exists(UPDATE)], collapse = ", "))
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)

# --- 1. Site cubes and the ones to patch ---------------------------------------
cubes <- list_site_cubes(ROOT)
wkt_files <- list.files(WKT_ROOT, pattern = "\\.wkt$", recursive = TRUE, full.names = TRUE)
cat(sprintf("%d site cube files, %d outline files.\n", nrow(cubes), length(wkt_files)))

# Misplaced cubes in the "unknown" folder: target path in the country's continent folder.
unk <- list.files(file.path(ROOT, "ramsar_site_data_100m_unknown"), pattern = "_data\\.csv$", recursive = TRUE, full.names = TRUE)
misplaced <- NULL
for (f in unk) {
  sid <- site_id_of(f); w <- find_site_wkt(sid, wkt_files)
  if (is.na(sid) || is.na(w)) { cat("  Cannot place", f, "- left as it is.\n"); next }
  country <- basename(dirname(w))
  conts <- CONTINENTS[dir.exists(file.path(ROOT, paste0("ramsar_site_data_100m_", CONTINENTS), country))]
  if (length(conts) != 1) { cat(sprintf("  %s: country %s is in %d continent folders - left as it is.\n", f, country, length(conts))); next }
  dest <- file.path(ROOT, paste0("ramsar_site_data_100m_", conts), country, paste0(sub("\\.wkt$", "", basename(w)), "_data.csv"))
  if (file.exists(dest)) { cat("  ", f, "already has a cube at", dest, "- left as it is.\n"); next }
  misplaced <- rbind(misplaced, data.frame(continent = conts, country = country, site_id = sid, file_path = dest,
                                           source_file = f, stringsAsFactors = FALSE))
}

if (has("--all")) {
  targets <- cubes
  targets$reason <- rep("all", nrow(targets))
} else if (!is.null(SITES)) {
  wanted <- if (file.exists(SITES)) trimws(readLines(SITES)) else trimws(strsplit(SITES, ",")[[1]])
  wanted <- wanted[nzchar(wanted)]
  targets <- cubes[cubes$site_id %in% wanted, ]
  targets$reason <- rep("requested", nrow(targets))
  if (length(setdiff(wanted, targets$site_id))) cat("Not found:", paste(setdiff(wanted, targets$site_id), collapse = " "), "\n")
} else {
  cat("Checking which cubes lack the update...\n")
  need <- vapply(cubes$file_path, needs_patch, logical(1), from_year = FROM)
  targets <- cubes[need, ]
  targets$reason <- rep("lacks update", nrow(targets))
}
targets$source_file <- targets$file_path
if (!is.null(misplaced) && is.null(SITES)) targets <- rbind(targets, cbind(misplaced, reason = "misplaced in unknown folder"))
targets$wkt <- if (nrow(targets)) unname(mapply(find_site_wkt, targets$site_id, country = targets$country,
                                                 MoreArgs = list(wkt_files = wkt_files))) else character()
if (any(is.na(targets$wkt))) {
  cat("No outline found, skipped:", paste(targets$site_id[is.na(targets$wkt)], collapse = " "), "\n")
  targets <- targets[!is.na(targets$wkt), ]
}
fwrite(targets, file.path(WORK, "patch_targets.csv"))
cat(sprintf("%d site cubes to patch:\n", nrow(targets)))
if (nrow(targets)) print(as.data.table(targets)[, .N, by = .(continent, country)][order(-N)][1:min(.N, 25)])
if (has("--dry-run")) quit(save = "no")
if (!nrow(targets)) { cat("Nothing to do.\n"); quit(save = "no") }

# --- 2. Outlines, 100 km squares -------------------------------------------------
outlines <- read_site_outlines(targets$site_id, targets$wkt)
sq <- unique(unlist(lapply(seq_len(nrow(outlines)), function(i) site_squares(sf::st_geometry(outlines)[i]))))
cat(sprintf("%d MGRS 100 km squares to search.\n", length(sq)))

# --- 3. Stream the update file(s) ---------------------------------------------------
# Rows from each update file go to their own folder (site_rows/1, site_rows/2, ...).
tmp_dir <- file.path(WORK, "site_rows"); dir.create(tmp_dir, showWarnings = FALSE)
old_tmp <- list.files(tmp_dir, recursive = TRUE, full.names = TRUE)
if (length(old_tmp)) stop("Work folder ", tmp_dir, " is not empty (from an earlier run?). Empty it and run again.")
cell_site <- data.table(mgrscellcode = character(), site_id = character())
seen <- data.table(mgrscellcode = character(), key = "mgrscellcode")   # cells already assigned
read_stats <- list()
for (k in seq_along(UPDATE)) {
  kdir <- file.path(tmp_dir, k); dir.create(kdir, showWarnings = FALSE)
  cat(sprintf("Reading update file %d of %d: %s (taxa: %s; years %d to %s)\n", k, length(UPDATE), UPDATE[k],
              paste(GROUPS[[k]], collapse = "+"), FROM, format(TO)))
  st <- c(rows_in_squares = 0, rows_other_years = 0, rows_other_taxa = 0, rows_kept = 0)
  n_lines <- stream_update(UPDATE[k], sq, chunk_lines = CHUNK, on_chunk = function(dt) {
    st["rows_in_squares"] <<- st["rows_in_squares"] + nrow(dt)
    y <- suppressWarnings(as.numeric(dt$year))
    in_years <- !is.na(y) & y >= FROM & y <= TO
    in_taxa <- if ("all" %in% GROUPS[[k]]) rep(TRUE, nrow(dt)) else in_subsets(dt, GROUPS[[k]])
    st["rows_other_years"] <<- st["rows_other_years"] + sum(!in_years)
    st["rows_other_taxa"] <<- st["rows_other_taxa"] + sum(in_years & !in_taxa)
    dt <- dt[in_years & in_taxa]
    st["rows_kept"] <<- st["rows_kept"] + nrow(dt)
    if (!nrow(dt)) return(invisible())
    cells <- data.table(mgrscellcode = unique(dt$mgrscellcode))
    new_cells <- cells[!seen, on = "mgrscellcode"]$mgrscellcode
    if (length(new_cells)) {
      cell_site <<- rbind(cell_site, assign_cells(new_cells, outlines))
      seen <<- rbind(seen, data.table(mgrscellcode = new_cells)); setkey(seen, mgrscellcode)
    }
    m <- merge(dt, cell_site, by = "mgrscellcode", allow.cartesian = TRUE)
    for (s in unique(m$site_id)) {
      f <- file.path(kdir, paste0(s, ".csv"))
      fwrite(m[m$site_id == s, CUBE_COLS, with = FALSE], f, append = file.exists(f))
    }
  })
  cat(sprintf("Read %s lines from %s. In the sites' 100 km squares: %s rows, of which %s from other years and %s from other taxa were left out.\n",
              format(n_lines, big.mark = ","), UPDATE[k], format(st["rows_in_squares"], big.mark = ","),
              format(st["rows_other_years"], big.mark = ","), format(st["rows_other_taxa"], big.mark = ",")))
  read_stats[[k]] <- data.table(file = UPDATE[k], taxa = paste(GROUPS[[k]], collapse = "+"), lines = n_lines, t(st))
}
fwrite(rbindlist(read_stats), file.path(WORK, "update_files_read.csv"))
empty_cube <- function() data.table(matrix(character(), 0, length(CUBE_COLS), dimnames = list(NULL, CUBE_COLS)))
read_new <- function(site_id, k) {
  nf <- file.path(tmp_dir, k, paste0(site_id, ".csv"))
  if (file.exists(nf)) fread(nf, colClasses = "character", showProgress = FALSE) else empty_cube()
}

# --- 4. Rebuild each site cube -------------------------------------------------------
summ <- list()
for (i in seq_len(nrow(targets))) {
  t <- targets[i, ]
  src <- t$source_file
  old <- if (file.exists(src)) fread(src, colClasses = "character", showProgress = FALSE) else data.table()
  if (t$reason == "misplaced in unknown folder") old <- old[0]       # update rows only; rebuilt below
  new_by_file <- lapply(seq_along(UPDATE), function(k) read_new(t$site_id, k))
  new <- rbindlist(new_by_file)
  if (!nrow(old)) old <- empty_cube()
  if (!nrow(new) && !sum(suppressWarnings(as.integer(old$year)) >= FROM, na.rm = TRUE) && t$reason != "misplaced in unknown folder") {
    summ[[i]] <- data.table(site_id = t$site_id, continent = t$continent, country = t$country, reason = "unchanged (no rows from the update year on, before or after)",
                            file = t$file_path, rows_before = nrow(old), rows_before_from = 0L, rows_update = 0L, rows_after = nrow(old),
                            rows_update_extra = 0L, occ_before_from = 0, occ_after_from = 0, reptile_rows_added = 0L, same_as_existing = NA)
    next
  }
  res <- patched_rows(old, new, FROM)
  # For cubes that already had the update: are the rows rebuilt from the first
  # update file the same as before, apart from reptiles (dropped by the June 2026
  # update; see the header)?
  had_update <- nrow(old) && !needs_patch_rows(old, FROM)
  main <- patched_rows(old, new_by_file[[1]], FROM)
  same_as_existing <- if (had_update) identical(row_keys(old[old$class != "Reptilia"], FROM, TO), row_keys(main[main$class != "Reptilia"], FROM, TO)) else NA
  dest <- if (is.null(OUT_DIR)) t$file_path else file.path(OUT_DIR, sub(paste0("^", ROOT, "/?"), "", t$file_path))
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  if (is.null(OUT_DIR) && file.exists(src)) {
    b <- file.path(BACKUP, sub(paste0("^", ROOT, "/?"), "", src))
    dir.create(dirname(b), showWarnings = FALSE, recursive = TRUE)
    if (!file.exists(b) && !file.copy(src, b, copy.date = TRUE)) stop("Backup failed for ", src)
  }
  tmpf <- paste0(dest, ".tmp")
  fwrite(res, tmpf)
  if (file.exists(dest)) file.remove(dest)
  if (!file.rename(tmpf, dest)) stop("Could not write ", dest)
  if (is.null(OUT_DIR) && t$reason == "misplaced in unknown folder") file.remove(src)   # a copy is in the backup folder
  yr <- function(d) { y <- suppressWarnings(as.integer(d$year)); y }
  summ[[i]] <- data.table(site_id = t$site_id, continent = t$continent, country = t$country, reason = t$reason,
                          file = dest, rows_before = nrow(old), rows_before_from = sum(yr(old) >= FROM, na.rm = TRUE),
                          rows_update = nrow(new), rows_update_extra = nrow(new) - nrow(new_by_file[[1]]), rows_after = nrow(res),
                          occ_before_from = sum(as.numeric(old$occurrences[yr(old) >= FROM]), na.rm = TRUE),
                          occ_after_from = sum(as.numeric(res$occurrences[yr(res) >= FROM]), na.rm = TRUE),
                          reptile_rows_added = sum(res$class == "Reptilia" & yr(res) >= FROM) - sum(old$class == "Reptilia" & yr(old) >= FROM),
                          same_as_existing = same_as_existing)
  if (i %% 25 == 0) cat(sprintf("  %d of %d sites written\n", i, nrow(targets)))
}
summ <- rbindlist(summ)
fwrite(summ, file.path(WORK, "patch_summary.csv"))
cat(sprintf("\nPatched %d site cubes%s. Rows from %d: %s before, %s after.\n", nrow(summ),
            if (is.null(OUT_DIR)) paste0(" in place (originals in ", BACKUP, ")") else paste0(" into ", OUT_DIR),
            FROM, format(sum(summ$rows_before_from), big.mark = ","), format(sum(summ$rows_update), big.mark = ",")))
if (any(!is.na(summ$same_as_existing))) {
  cat(sprintf("Cubes that already had the update: %d; rebuilt rows identical apart from reptiles: %d; reptile rows added: %s.\n",
              sum(!is.na(summ$same_as_existing)), sum(summ$same_as_existing, na.rm = TRUE),
              format(sum(summ$reptile_rows_added[!is.na(summ$same_as_existing)]), big.mark = ",")))
  if (any(summ$same_as_existing == FALSE, na.rm = TRUE)) print(summ[same_as_existing == FALSE, .(site_id, rows_before_from, rows_update)])
}
cat("Summary:", file.path(WORK, "patch_summary.csv"), "\n")
if (length(UPDATE) > 1) cat(sprintf("Rows from the second and later update files: %s.\n", format(sum(summ$rows_update_extra), big.mark = ",")))
cat("The per-site rows in", tmp_dir, "can be deleted once the result is checked.\n")
