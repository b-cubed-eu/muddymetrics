#!/usr/bin/env Rscript
# Checks every site cube against its sources before any analysis is run.
#
# From the source files alone, this rebuilds what each site cube should
# contain and compares it row by row with the cube on disk:
#
#   - outlines: the Ramsar shapefile (inst/extdata/Ramsar_boundaries_160626),
#     not the .wkt files the cubes were split with; the .wkt files are
#     compared with the shapefile separately;
#   - years before --from-year (2024): the continental cubes
#     (inst/extdata/continental_gbif_data/*_100m_mgrs_all.csv), all of them
#     for every site, so data filed under a neighbouring continent is found too;
#   - years from --from-year on: the update files given with --update (the
#     2024-2026 update and the fish download).
#
# A cell belongs to a site when the point mgrs::mgrs_to_latlng() returns for it
# (the south-west corner of the 100 m cell) lies in or on the site outline
# (planar geometry), the rule the cubes were built with. Each check is
# explained in the report the last step writes.
#
# The code is separate from the code that built the cubes (see
# scripts/lib/cube_check.R), and the script never deletes or changes a cube.
#
# Usage (from the repository root):
#   Rscript scripts/check_site_cubes.R --step all --update <main update>,<fish download>
#   Rscript scripts/check_site_cubes.R --step inventory          # files and outlines only (minutes)
#   Rscript scripts/check_site_cubes.R --step reference --update <main>,<fish> --cores 8
#        # reads all sources once (hours); resumable, each finished byte range is kept
#   Rscript scripts/check_site_cubes.R --step compare --tag after_patch
#   Rscript scripts/check_site_cubes.R --step compare --tag after_taxonomy --name-table _local_data/taxonomy/name_table.csv
#        # after harmonise_cube_taxonomy.R --apply: the reference is renamed with the name table first
#   Rscript scripts/check_site_cubes.R --step report --tag after_patch
#
# Options:
#   --step          inventory | reference | compare | report | all
#   --out           output folder [_local_data/cube_check]
#   --root          data folder [inst/extdata]
#   --shapefile     [inst/extdata/Ramsar_boundaries_160626/features_publishedPolygon.shp]
#   --wkt-dir       [inst/extdata/ramsar_sites_wkt]
#   --continental   folder with the continental cubes [inst/extdata/continental_gbif_data]
#   --update        comma-separated update files [the 2024-2026 update in
#                   inst/extdata/ramsar_site_data_100m_world_2024_2026]
#   --from-year     first year taken from the update files [2024]
#   --cores         parallel processes for the reference step [6]
#   --compare-cores parallel processes for the compare step [4]
#   --part-gb       size of the byte ranges sources are read in [2]
#   --chunk-lines   lines read at a time [1000000]
#   --buckets       number of files the reference rows are spread over [200]
#   --tag           name of the compare run (output in <out>/compare_<tag>) [check]; buckets
#                   whose cubes changed since the last run with the same tag are compared again
#   --name-table    compare: harmonised cubes; apply this name table to the reference
#   --sites         compare: only these site ids (comma-separated or a file)
#
# Output in --out: inventory.csv, outline_check.csv, outlines.rds,
# reference/ (rows per byte range and bucket, sources.csv), compare_<tag>/
# (sites.csv, year_group.csv.gz, examples.csv) and REPORT_<tag>.md.

suppressPackageStartupMessages({ library(data.table); library(sf); library(mgrs); library(stringi); library(parallel) })
source("scripts/lib/sufficiency_criteria.R")   # SUBSETS and subset_mask(), for the taxon groups in the report
source("scripts/lib/cube_check.R")

arg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
STEP <- match.arg(arg("--step", "all"), c("inventory", "reference", "compare", "report", "all"))
OUT <- arg("--out", "_local_data/cube_check")
ROOT <- arg("--root", "inst/extdata")
SHP <- arg("--shapefile", file.path(ROOT, "Ramsar_boundaries_160626", "features_publishedPolygon.shp"))
WKT_DIR <- arg("--wkt-dir", file.path(ROOT, "ramsar_sites_wkt"))
CONT_DIR <- arg("--continental", file.path(ROOT, "continental_gbif_data"))
UPDATE <- trimws(strsplit(arg("--update", file.path(ROOT, "ramsar_site_data_100m_world_2024_2026", "0074420-260519110011954.csv")), ",")[[1]])
FROM <- as.integer(arg("--from-year", "2024"))
CORES <- as.integer(arg("--cores", "6")); CCORES <- as.integer(arg("--compare-cores", "4"))
PART_BYTES <- as.numeric(arg("--part-gb", "2")) * 1024^3
CHUNK <- as.numeric(arg("--chunk-lines", "1e6"))
NB <- as.integer(arg("--buckets", "200"))
TAG <- arg("--tag", "check")
NAME_TABLE <- arg("--name-table")
SITES <- arg("--sites")
REF_DIR <- file.path(OUT, "reference"); CMP_DIR <- file.path(OUT, paste0("compare_", TAG))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
settings <- data.table(setting = c("out", "shapefile", "wkt_dir", "continental", "update", "from_year", "buckets", "tag", "name_table"),
                       value = c(OUT, SHP, WKT_DIR, CONT_DIR, paste(UPDATE, collapse = ","), FROM, NB, TAG, if (is.null(NAME_TABLE)) "" else NAME_TABLE))
say <- function(...) cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")

# --- inventory: files, outlines ------------------------------------------------
if (STEP %in% c("inventory", "all")) {
  say("Reading the Ramsar shapefile %s", SHP)
  outlines <- chk_read_outlines(SHP)
  saveRDS(outlines, file.path(OUT, "outlines.rds"))
  cubes <- chk_list_cubes(ROOT)
  wkt <- list.files(WKT_DIR, pattern = "\\.wkt$", recursive = TRUE, full.names = TRUE)
  say("%d sites in the shapefile, %d cube files, %d outline files", nrow(outlines), nrow(cubes), length(wkt))
  inv <- chk_inventory(outlines, cubes, wkt)
  fwrite(inv, file.path(OUT, "inventory.csv"))
  fwrite(cubes, file.path(OUT, "cube_files.csv"))
  say("Comparing the outline files with the shapefile...")
  oc <- chk_compare_wkt(outlines, wkt)
  fwrite(oc, file.path(OUT, "outline_check.csv"))
  fwrite(settings, file.path(OUT, "settings_inventory.csv"))
  say("Inventory: %d sites with an issue; outline files more than 10 m from the shapefile on average: %d; not compared: %d",
      sum(inv$issues != ""), sum(oc$mean_shift_m > 10, na.rm = TRUE), sum(oc$status != "compared"))
}

# --- reference: read all sources -----------------------------------------------
if (STEP %in% c("reference", "all")) {
  if (!file.exists(file.path(OUT, "outlines.rds"))) stop("Run --step inventory first.")
  cont <- list.files(CONT_DIR, pattern = "_100m_mgrs_all\\.csv$", full.names = TRUE)
  if (!length(cont)) stop("No continental cubes (*_100m_mgrs_all.csv) in ", CONT_DIR)
  if (any(!file.exists(UPDATE))) stop("Update file not found: ", paste(UPDATE[!file.exists(UPDATE)], collapse = ", "))
  src <- rbind(data.table(path = cont, role = "before"), data.table(path = UPDATE, role = "from"))
  plan <- rbindlist(lapply(seq_len(nrow(src)), function(i) cbind(chk_plan_parts(src$path[i], PART_BYTES), role = src$role[i])))
  plan[, part_dir := file.path(REF_DIR, sprintf("%s_part%04d", tools::file_path_sans_ext(basename(path)), part))]
  plan[, eol_bytes := vapply(path, chk_eol_bytes, 1L)]
  plan[, final_eol := vapply(path, chk_final_eol, TRUE)]
  dir.create(REF_DIR, showWarnings = FALSE, recursive = TRUE)
  todo <- plan[!file.exists(file.path(part_dir, "done.csv"))]
  say("Sources: %d files, %.1f GB, %d byte ranges (%d still to read) on %d cores",
      nrow(src), sum(file.size(src$path)) / 1024^3, nrow(plan), nrow(todo), CORES)
  if (nrow(todo)) {
    cl <- makePSOCKcluster(CORES)
    res <- tryCatch({
      clusterExport(cl, c("OUT", "FROM", "NB", "CHUNK"))
      clusterEvalQ(cl, {
        suppressPackageStartupMessages({ library(data.table); library(sf); library(mgrs); library(stringi) })
        source("scripts/lib/cube_check.R")
        outlines <- readRDS(file.path(OUT, "outlines.rds"))
        NULL
      })
      parLapplyLB(cl, seq_len(nrow(todo)), function(i, todo) {
        t <- todo[i]
        tryCatch(chk_reference_part(t$path, t$role, t$part, t$start, t$end, outlines, FROM, t$part_dir, NB, CHUNK),
                 error = function(e) data.table(source = t$path, part = t$part, error = conditionMessage(e)))
      }, todo = todo, chunk.size = 1)
    }, finally = stopCluster(cl))
    err <- rbindlist(res, fill = TRUE)
    if ("error" %in% names(err) && any(!is.na(err$error))) {
      print(err[!is.na(error), .(source, part, error)])
      stop("Some byte ranges failed (listed above). Fix the cause and run the step again; finished ranges are kept.")
    }
  }
  parts <- rbindlist(lapply(plan$part_dir, function(d) fread(file.path(d, "done.csv"))), fill = TRUE)
  parts <- merge(parts, unique(plan[, .(source = path, part, size, eol_bytes, final_eol)]), by = c("source", "part"))
  til <- parts[, { t <- chk_check_tiling(.SD); .(read_completely = t$ok, gaps_after_part = paste(t$gaps_after_part, collapse = ","),
                                               bytes_read_to = t$end_read, file_end = t$end_expected) }, by = source]
  fwrite(parts, file.path(REF_DIR, "sources_parts.csv"))
  fwrite(merge(parts[, lapply(.SD, sum), by = .(source, role),
                     .SDcols = c("lines", "n_empty", "n_wrong_fields", "n_bad_utf8", "n_invalid_code", "n_no_coords",
                                 "rows_in_sites", "rows_kept", "rows_other_period", "minutes")], til, by = "source"),
         file.path(REF_DIR, "sources.csv"))
  fwrite(settings, file.path(REF_DIR, "settings.csv"))
  say("All sources read; read completely: %s", paste(til$read_completely, collapse = ","))
  if (!all(til$read_completely)) warning("A source file was not read completely: see ", file.path(REF_DIR, "sources.csv"))
}

# --- compare: cubes against the reference ----------------------------------------
if (STEP %in% c("compare", "all")) {
  if (!file.exists(file.path(REF_DIR, "sources.csv"))) stop("Run --step reference first.")
  dir.create(CMP_DIR, showWarnings = FALSE, recursive = TRUE)
  outlines <- readRDS(file.path(OUT, "outlines.rds"))
  cubes <- chk_list_cubes(ROOT)
  ids <- unique(c(outlines$site_id, cubes$site_id[!is.na(cubes$site_id)]))
  if (!is.null(SITES)) ids <- intersect(ids, if (file.exists(SITES)) trimws(readLines(SITES)) else trimws(strsplit(SITES, ",")[[1]]))
  rs <- fread(file.path(REF_DIR, "sources_parts.csv"))
  if (any(rs$n_buckets != NB) || any(rs$from_year != FROM))
    stop("The reference was built with --buckets ", paste(unique(rs$n_buckets), collapse = ","), " and --from-year ",
         paste(unique(rs$from_year), collapse = ","), "; use the same settings for the compare step.")
  buckets <- sort(unique(chk_bucket_of(ids, NB)))
  # A bucket is redone when it has no result yet, when its result lacks some of the sites asked for
  # (an earlier run with --sites), or when one of its cubes changed after the result was written.
  redo <- vapply(buckets, function(b) {
    f <- file.path(CMP_DIR, sprintf("bucket_%03d_sites.csv", b))
    if (!file.exists(f)) return(TRUE)
    want <- ids[chk_bucket_of(ids, NB) == b]
    if (!all(want %in% fread(f, select = "site_id", colClasses = "character")$site_id)) return(TRUE)
    cf <- cubes$file[!is.na(cubes$site_id) & cubes$site_id %in% want]
    length(cf) > 0 && max(as.numeric(file.mtime(cf))) >= as.numeric(file.mtime(f))
  }, TRUE)
  todo <- buckets[redo]
  say("Comparing %d sites in %d buckets (%d to do) on %d cores%s", length(ids), length(buckets), length(todo), CCORES,
      if (is.null(NAME_TABLE)) "" else paste0(", reference renamed with ", NAME_TABLE))
  run_bucket <- function(b) {
    nt <- if (!is.null(NAME_TABLE)) fread(NAME_TABLE, colClasses = "character") else NULL
    cols <- if (is.null(nt)) CHK_COLS else CHK_HARMONISED_COLS
    ref_all <- chk_read_reference_bucket(REF_DIR, b)
    sites <- ids[chk_bucket_of(ids, NB) == b]
    res <- lapply(sites, function(s) {
      files <- cubes$file[cubes$site_id %in% s]
      cube <- if (length(files)) chk_read_cube(files) else data.table(matrix(character(), 0, length(CHK_COLS), dimnames = list(NULL, CHK_COLS)))
      cube <- chk_normalise(cube, CHK_COLS)
      ref <- chk_normalise(ref_all[site_id == s], CHK_COLS)
      internal <- chk_cube_internal(cube, FROM)
      yg <- chk_year_group_table(cube, ref, SUBSETS)
      if (!is.null(nt)) { ref <- chk_apply_name_table(ref[, CHK_HARMONISED_COLS, with = FALSE], nt); cube <- cube[, CHK_HARMONISED_COLS, with = FALSE] }
      j <- match(s, outlines$site_id)
      cmp <- if (is.na(j)) list(counts = data.table(rows_cube = nrow(cube), rows_ref = nrow(ref)), examples = NULL) else
        chk_compare_site(cube, ref, cols, outlines[j, ], FROM)
      list(site = cbind(data.table(site_id = s, in_shapefile = !is.na(j), n_cube_files = length(files)), cmp$counts, internal),
           yg = if (nrow(yg)) cbind(site_id = s, yg) else NULL,
           ex = if (!is.null(cmp$examples) && nrow(cmp$examples)) cbind(site_id = s, cmp$examples) else NULL)
    })
    fwrite(rbindlist(lapply(res, `[[`, "yg")), file.path(CMP_DIR, sprintf("bucket_%03d_year_group.csv", b)))
    fwrite(rbindlist(lapply(res, `[[`, "ex"), fill = TRUE), file.path(CMP_DIR, sprintf("bucket_%03d_examples.csv", b)))
    f <- file.path(CMP_DIR, sprintf("bucket_%03d_sites.csv", b))   # written last: marks the bucket as done
    fwrite(rbindlist(lapply(res, `[[`, "site"), fill = TRUE), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f)
    b
  }
  if (length(todo)) {
    cl <- makePSOCKcluster(CCORES)
    res <- tryCatch({
      clusterExport(cl, c("OUT", "REF_DIR", "CMP_DIR", "FROM", "NB", "NAME_TABLE", "ids", "cubes", "run_bucket"), envir = environment())
      clusterEvalQ(cl, {
        suppressPackageStartupMessages({ library(data.table); library(sf); library(mgrs); library(stringi) })
        source("scripts/lib/sufficiency_criteria.R"); source("scripts/lib/cube_check.R")
        outlines <- readRDS(file.path(OUT, "outlines.rds"))
        NULL
      })
      parLapplyLB(cl, todo, function(b) tryCatch(run_bucket(b), error = function(e) paste("bucket", b, ":", conditionMessage(e))), chunk.size = 1)
    }, finally = stopCluster(cl))
    errs <- unlist(res[vapply(res, is.character, TRUE)])
    if (length(errs)) { print(errs); stop("Some buckets failed (listed above); finished buckets are kept.") }
  }
  rd <- function(pat) rbindlist(lapply(list.files(CMP_DIR, pattern = pat, full.names = TRUE), fread, colClasses = list(character = "site_id")), fill = TRUE)
  sites <- rd("^bucket_[0-9]+_sites\\.csv$")
  fwrite(sites[order(as.integer(sub("site_", "", site_id)))], file.path(CMP_DIR, "sites.csv"))
  fwrite(rd("^bucket_[0-9]+_year_group\\.csv$"), file.path(CMP_DIR, "year_group.csv.gz"))
  fwrite(rd("^bucket_[0-9]+_examples\\.csv$"), file.path(CMP_DIR, "examples.csv"))
  fwrite(settings, file.path(CMP_DIR, "settings.csv"))
  say("Compared %d sites: %d identical to their sources", nrow(sites),
      sum(sites$rows_cube == sites$rows_same & sites$rows_ref == sites$rows_same, na.rm = TRUE))
}

# --- report ----------------------------------------------------------------------
if (STEP %in% c("report", "all")) {
  source("scripts/lib/cube_check_report.R")
  f <- chk_write_report(OUT, TAG, FROM)
  say("Report written: %s", f)
}
