# Matching old-backbone species names in the site cubes to the new GBIF backbone.
#
# Rows up to 2023 in the site cubes carry the old GBIF backbone (numeric keys);
# rows from the June 2026 update carry the new backbone, which is built from the
# Catalogue of Life (checklist 7ddf754f-d193-4cc9-b351-99906754a03b; short
# alphanumeric keys such as "7CYVX"). Some species have different names in the
# two (Sylvia curruca -> Curruca curruca), so counting species by name across
# years treats them as two species.
#
# The name table maps each old (kingdom, species) to the accepted species of the
# new backbone with its full classification:
#   1. names that occur as species in the new-backbone rows are taken from there;
#   2. all others are matched with the GBIF species-match API (v2) against the
#      new backbone. Synonyms resolve to their accepted name.
# Used by scripts/harmonise_cube_taxonomy.R and its tests.

NEW_BACKBONE <- "7ddf754f-d193-4cc9-b351-99906754a03b"
TAXON_COLS <- c("kingdom", "kingdomkey", "phylum", "phylumkey", "class", "classkey",
                "order", "orderkey", "family", "familykey", "genus", "genuskey",
                "species", "specieskey")
RANKS <- c("KINGDOM", "PHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES")

is_old_key <- function(k) grepl("^[0-9]+$", k)

#' URL for one name in the GBIF v2 species-match API.
match_url <- function(species, kingdom = NA, phylum = NA, class = NA, checklist = NEW_BACKBONE) {
  one <- function(s, k, p, c) {
    q <- c(checklistKey = checklist, scientificName = s, kingdom = k, phylum = p, class = c)
    q <- q[!is.na(q) & nzchar(q)]
    paste0("https://api.gbif.org/v2/species/match?",
           paste(names(q), vapply(q, utils::URLencode, "", reserved = TRUE), sep = "=", collapse = "&"))
  }
  n <- length(species)
  unname(mapply(one, species, rep_len(kingdom, n), rep_len(phylum, n), rep_len(class, n)))
}

#' Turn one match response (JSON text) into a one-row data.frame with the
#' accepted species and its classification in cube column names, plus the
#' match type, confidence and whether the old name is a synonym.
#' Rows that do not resolve to a species get status "unmatched".
parse_match <- function(json) {
  empty <- as.list(setNames(rep(NA_character_, length(TAXON_COLS)), TAXON_COLS))
  out <- function(status, match_type = NA, confidence = NA, synonym = NA, cls = empty)
    data.frame(c(cls, list(status = status, match_type = match_type, confidence = as.numeric(confidence),
                           synonym = synonym)), stringsAsFactors = FALSE)
  j <- tryCatch(jsonlite::fromJSON(json, simplifyVector = FALSE), error = function(e) NULL)
  if (is.null(j) || is.null(j$usage)) return(out("unmatched"))
  mt <- j$diagnostics$matchType %||% NA; cf <- j$diagnostics$confidence %||% NA
  syn <- isTRUE(j$synonym)
  cl <- j$classification
  if (is.null(cl) || !length(cl)) return(out("unmatched", mt, cf, syn))
  ranks <- vapply(cl, function(x) x$rank %||% "", ""); names_ <- vapply(cl, function(x) x$name %||% NA_character_, "")
  keys <- vapply(cl, function(x) x$key %||% NA_character_, "")
  if (!"SPECIES" %in% ranks || !mt %in% c("EXACT", "FUZZY", "VARIANT")) return(out("unmatched", mt, cf, syn))
  cls <- empty
  for (r in RANKS) {
    i <- which(ranks == r)[1]
    if (is.na(i)) next
    nm <- tolower(r); cls[[nm]] <- names_[i]; cls[[paste0(nm, "key")]] <- keys[i]
  }
  out(if (syn) "synonym" else "accepted", mt, cf, syn, cls)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Fetch many URLs concurrently with the curl package; returns the response
#' bodies as text (NA where a request failed after `retries` attempts).
fetch_urls <- function(urls, connections = 4, retries = 3, timeout = 60) {
  res <- rep(NA_character_, length(urls))
  todo <- seq_along(urls)
  for (attempt in seq_len(retries)) {
    if (!length(todo)) break
    pool <- curl::new_pool(total_con = connections, host_con = connections)
    for (i in todo) local({
      ii <- i
      h <- curl::new_handle(timeout = timeout, useragent = "muddymetrics taxonomy harmonisation (R curl)")
      curl::curl_fetch_multi(urls[ii], handle = h, pool = pool,
        done = function(r) if (r$status_code == 200) res[ii] <<- rawToChar(r$content),
        fail = function(msg) NULL)
    })
    curl::multi_run(pool = pool)
    todo <- which(is.na(res))
    if (length(todo)) Sys.sleep(5 * attempt)
  }
  res
}

#' Name table for old (kingdom, species) pairs.
#' `old`: data.frame with kingdom, phylum, class, species (one row per name).
#' `new_rows`: data.frame of new-backbone taxa (TAXON_COLS) seen in the cubes.
#' `cache`: path of a CSV holding earlier API results; new results are appended,
#' so an interrupted run can be continued. `fetch`: function(urls) -> bodies.
#' `limit`: only query this many names via the API (for a first test).
build_name_table <- function(old, new_rows, cache, fetch = fetch_urls, batch = 500, limit = Inf, progress = TRUE) {
  old <- unique(data.table::as.data.table(old)[, c("kingdom", "phylum", "class", "species"), with = FALSE], by = c("kingdom", "species"))
  nr <- unique(data.table::as.data.table(new_rows)[, TAXON_COLS, with = FALSE], by = c("kingdom", "species"))
  direct <- merge(old[, .(kingdom, species)], nr, by = c("kingdom", "species"))
  direct[, `:=`(old_species = species, old_kingdom = kingdom, status = "in update data", match_type = NA_character_,
                confidence = NA_real_, synonym = NA)]
  rest <- old[!direct, on = c("kingdom", "species")]
  done <- if (file.exists(cache)) data.table::fread(cache, colClasses = "character") else NULL
  if (!is.null(done) && nrow(done)) rest <- rest[!done, on = c(kingdom = "old_kingdom", species = "old_species")]
  if (is.finite(limit)) rest <- rest[seq_len(min(nrow(rest), limit))]
  if (progress) cat(sprintf("%d names from the update data, %d already matched earlier, %d to match via the API.\n",
                            nrow(direct), if (is.null(done)) 0L else nrow(done), nrow(rest)))
  for (b in split(seq_len(nrow(rest)), ceiling(seq_len(nrow(rest)) / batch))) {
    r <- rest[b]
    body <- fetch(match_url(r$species, r$kingdom, r$phylum, r$class))
    p <- data.table::rbindlist(lapply(body, function(x) if (is.na(x)) parse_match("{}") else parse_match(x)))
    p[is.na(body), status := "request failed"]
    p[, `:=`(old_kingdom = r$kingdom, old_species = r$species)]
    p <- p[status != "request failed"]                       # retried on the next run
    data.table::fwrite(p, cache, append = file.exists(cache))
    if (progress) cat(sprintf("  %d of %d names matched via the API\n", max(b), nrow(rest)))
  }
  api <- if (file.exists(cache)) data.table::fread(cache, colClasses = "character") else NULL
  cols <- c("old_kingdom", "old_species", TAXON_COLS, "status", "match_type", "confidence", "synonym")
  rbind(direct[, cols, with = FALSE], if (!is.null(api)) api[, cols, with = FALSE], fill = TRUE)
}

#' Apply a name table to one cube: rows with old-backbone keys whose
#' (kingdom, species) resolved to a species get the new classification; rows
#' that then share cell, year and species key are merged (occurrences summed,
#' uncertainties take the minimum, other columns from the first row).
#' Returns list(rows, n_renamed, n_merged).
harmonise_rows <- function(dt, name_table) {
  dt <- data.table::as.data.table(dt)
  for (cc in names(dt)) data.table::set(dt, j = cc, value = as.character(dt[[cc]]))
  nt <- data.table::as.data.table(name_table)[status %in% c("in update data", "accepted", "synonym")]
  nt <- unique(nt, by = c("old_kingdom", "old_species"))
  dt[, .row := .I]
  old <- dt[is_old_key(specieskey)]
  m <- merge(old[, .(.row, old_kingdom = kingdom, old_species = species)], nt, by = c("old_kingdom", "old_species"))
  renamed <- sum(m$species != m$old_species, na.rm = TRUE)
  for (cc in TAXON_COLS) data.table::set(dt, i = m$.row, j = cc, value = m[[cc]])
  dt[, .row := NULL]
  key <- c("mgrscellcode", "year", "specieskey")
  dup <- duplicated(dt, by = key) | duplicated(dt, by = key, fromLast = TRUE)
  n_before <- nrow(dt)
  if (any(dup)) {
    num_min <- intersect(c("mintemporaluncertainty", "mincoordinateuncertaintyinmeters"), names(dt))
    d <- dt[dup]
    first <- d[!duplicated(d, by = key)]
    agg <- d[, c(list(occurrences = as.character(sum(as.numeric(occurrences), na.rm = TRUE))),
                 lapply(.SD, function(x) { v <- suppressWarnings(as.numeric(x)); if (all(is.na(v))) NA_character_ else as.character(min(v, na.rm = TRUE)) })),
             by = key, .SDcols = num_min]
    first[, c("occurrences", num_min) := NULL]
    first <- merge(first, agg, by = key, sort = FALSE)[, names(dt), with = FALSE]
    dt <- rbind(dt[!dup], first)
  }
  list(rows = dt, n_renamed = renamed, n_merged = n_before - nrow(dt))
}
