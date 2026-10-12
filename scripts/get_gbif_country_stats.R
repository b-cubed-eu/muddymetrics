# Country-level GBIF statistics for the model of which Ramsar Sites can be
# verified (Paper 1): how much georeferenced data GBIF holds for each
# Contracting Party, how much of it is published from within the country, how
# many publishers contribute, and the country's GBIF participation status.
#
# Usage (from the repository folder):
#   Rscript scripts/get_gbif_country_stats.R [--countries _local_data/RSIS/rsis_country_iso2.csv]
#           [--out _local_data/gbif_country] [--years 1990,2024] [--delay 0.2]
#           [--api https://api.gbif.org/v1] [--facet-limit 20000]
# Needs the jsonlite package. Every API response is kept in <out>/cache/ (one
# JSON file per request), so a rerun only requests what is missing; delete the
# cache to refresh. Writes <out>/gbif_country_stats.csv and <out>/gbif_nodes.csv.
#
# Per country (ISO 3166-1 alpha-2; territories such as Greenland or the UK
# overseas territories are not counted under their Party), for occurrences with
# coordinates and without geospatial issues in the years given:
#   occ_in_country           occurrences located in the country
#   occ_in_country_domestic  of these, published by an organisation in the country
#   occ_published            occurrences published from the country (anywhere)
#   publishers_in_country    publishing organisations with occurrences located in the country
#   publishers_domestic      organisations in the country that publish occurrences (anywhere)
#   orgs_registered          organisations registered in the country (all, and endorsed only)
#   participation_status     from the GBIF participant (node) list, country nodes only, e.g. VOTING, ASSOCIATE
# The publisher counts come from facets; *_capped is TRUE if a count reached
# the facet limit (then it is a lower bound).

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Please install the jsonlite package")
args <- commandArgs(TRUE)
opt <- function(name, default) { i <- match(paste0("--", name), args); if (is.na(i)) default else args[i + 1] }
COUNTRIES <- opt("countries", "_local_data/RSIS/rsis_country_iso2.csv")
OUT   <- opt("out", "_local_data/gbif_country")
YEARS <- opt("years", "1990,2024")
DELAY <- as.numeric(opt("delay", "0.2"))
API   <- sub("/$", "", opt("api", "https://api.gbif.org/v1"))
FACET_LIMIT <- as.integer(opt("facet-limit", "20000"))
CACHE <- file.path(OUT, "cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

# One GET request, cached on disk under a name derived from the URL.
get_json <- function(url) {
  f <- file.path(CACHE, paste0(gsub("[^A-Za-z0-9]+", "_", sub("^https?://", "", url)), ".json"))
  if (nchar(basename(f)) > 200) f <- file.path(CACHE, paste0(substr(basename(f), 1, 150), "_", nchar(url), ".json"))
  if (!file.exists(f)) {
    ok <- FALSE
    for (attempt in 1:3) {
      res <- tryCatch({ utils::download.file(url, paste0(f, ".part"), quiet = TRUE, mode = "wb", method = "libcurl"); TRUE },
                      error = function(e) { message("  retry ", attempt, ": ", conditionMessage(e)); FALSE })
      Sys.sleep(DELAY * attempt)
      if (res) { file.rename(paste0(f, ".part"), f); ok <- TRUE; break }
    }
    if (!ok) stop("Request failed three times: ", url)
  }
  jsonlite::fromJSON(f, simplifyVector = FALSE)
}
# For the optional counts: NA (with a message) instead of stopping.
safe <- function(expr) tryCatch(expr, error = function(e) { message("  ", conditionMessage(e)); NA })
occ_filter <- sprintf("hasCoordinate=true&hasGeospatialIssue=false&year=%s", YEARS)
occ_count <- function(q) get_json(sprintf("%s/occurrence/search?%s&%s&limit=0", API, q, occ_filter))$count
facet_n <- function(q, field) {
  j <- get_json(sprintf("%s/occurrence/search?%s&%s&limit=0&facet=%s&facetLimit=%d&facetMincount=1",
                        API, q, occ_filter, field, FACET_LIMIT))
  n <- if (length(j$facets)) length(j$facets[[1]]$counts) else 0L
  c(n = n, capped = n >= FACET_LIMIT)
}

# na.strings = "": Namibia's code is "NA".
cn <- utils::read.csv(COUNTRIES, encoding = "UTF-8", stringsAsFactors = FALSE, na.strings = "")
stopifnot(all(c("country", "iso2") %in% names(cn)), !anyDuplicated(cn$iso2), all(grepl("^[A-Z]{2}$", cn$iso2)))

# Participant (node) list, paged.
nodes <- list(); off <- 0
repeat {
  j <- get_json(sprintf("%s/node?limit=200&offset=%d", API, off))
  nodes <- c(nodes, j$results)
  if (isTRUE(j$endOfRecords) || !length(j$results)) break
  off <- off + 200
}
fld <- function(x, k) if (is.null(x[[k]])) NA_character_ else as.character(x[[k]])
nd <- data.frame(node_key = vapply(nodes, fld, "", "key"), title = vapply(nodes, fld, "", "title"),
                 type = vapply(nodes, fld, "", "type"), participation_status = vapply(nodes, fld, "", "participationStatus"),
                 country = vapply(nodes, fld, "", "country"), participant_since = vapply(nodes, fld, "", "participantSince"),
                 gbif_region = vapply(nodes, fld, "", "gbifRegion"), stringsAsFactors = FALSE)
utils::write.csv(nd, file.path(OUT, "gbif_nodes.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Nodes: %d (with a country: %d)\n", nrow(nd), sum(!is.na(nd$country))))

rows <- vector("list", nrow(cn))
for (i in seq_len(nrow(cn))) {
  cc <- cn$iso2[i]
  pin <- safe(facet_n(sprintf("country=%s", cc), "publishingOrg")); if (anyNA(pin)) pin <- c(n = NA, capped = NA)
  pdo <- safe(facet_n(sprintf("publishingCountry=%s", cc), "publishingOrg")); if (anyNA(pdo)) pdo <- c(n = NA, capped = NA)
  # Country nodes only: other node types (e.g. regional or thematic networks) also carry a country.
  nd_c <- nd[!is.na(nd$country) & nd$country == cc & nd$type %in% "COUNTRY", ]
  rows[[i]] <- data.frame(
    country = cn$country[i], iso2 = cc,
    occ_in_country = occ_count(sprintf("country=%s", cc)),
    occ_in_country_domestic = occ_count(sprintf("country=%s&publishingCountry=%s", cc, cc)),
    occ_published = occ_count(sprintf("publishingCountry=%s", cc)),
    publishers_in_country = pin[["n"]], publishers_in_country_capped = as.logical(pin[["capped"]]),
    publishers_domestic = pdo[["n"]], publishers_domestic_capped = as.logical(pdo[["capped"]]),
    orgs_registered = safe(get_json(sprintf("%s/organization?country=%s&limit=0", API, cc))$count),
    orgs_registered_endorsed = safe(get_json(sprintf("%s/organization?country=%s&isEndorsed=true&limit=0", API, cc))$count),
    participation_status = if (nrow(nd_c)) paste(sort(unique(nd_c$participation_status)), collapse = ";") else NA_character_,
    participant_since = if (nrow(nd_c)) paste(sort(unique(nd_c$participant_since)), collapse = ";") else NA_character_,
    stringsAsFactors = FALSE)
  if (i %% 20 == 0) cat(sprintf("  %d of %d countries\n", i, nrow(cn)))
}
st <- do.call(rbind, rows)
st$share_domestic <- ifelse(st$occ_in_country > 0, st$occ_in_country_domestic / st$occ_in_country, NA_real_)
st$years <- YEARS
st$retrieved <- format(Sys.Date())
utils::write.csv(st, file.path(OUT, "gbif_country_stats.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat(sprintf("\nCountries: %d; with GBIF participation status: %d; any capped publisher count: %s\n",
            nrow(st), sum(!is.na(st$participation_status)), any(st$publishers_in_country_capped | st$publishers_domestic_capped, na.rm = TRUE)))
cat(sprintf("Occurrences located in these countries (%s): %.0f million\n", YEARS, sum(st$occ_in_country) / 1e6))
cat(sprintf("Written: %s\n", file.path(OUT, "gbif_country_stats.csv")))
