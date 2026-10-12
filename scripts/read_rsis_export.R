# Read the tabular export of the Ramsar Sites Information Service (RSIS,
# "ramsar-sites-tabular-information.csv", one row per Ramsar Site) into two
# tidy tables for the analysis.
#
# Usage: Rscript read_rsis_export.R <export.csv> <out_dir>
# Writes:
#   rsis_sites.csv          one row per site (see the column list at the end)
#   rsis_wetland_types.csv  one row per site and wetland type: group (marine,
#                           inland, human_made), Ramsar code, label
# Stops if the export is incomplete (rows with the wrong number of fields,
# duplicated or non-numeric site numbers, unparseable dates or areas).
suppressPackageStartupMessages(library(data.table))
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
args <- commandArgs(TRUE)
x <- fread(args[1], encoding = "UTF-8", colClasses = "character", fill = FALSE)
out_dir <- args[2]

# --- checks on the export ------------------------------------------------------
expected <- c("Ramsar Site No.", "Site name", "Region", "Country", "Territory", "Designation date",
              "Last publication date", "Area (ha)", "Latitude", "Longitude", "Annotated summary",
              paste0("Criterion", 1:9), "Wetland Type", "Maximum elevation", "Minimum elevation",
              "Montreux listed", "Management plan implemented", "Management plan available",
              "Ecosystem services", "Threats", "large administrative region",
              "Global international legal designations", "Regional international legal designations",
              "National conservation designation",
              "Does the wetland extend onto the territory of one or more other countries?",
              "Ramsar Advisory Mission?")
stopifnot(identical(names(x), expected))
stopifnot(all(grepl("^[0-9]+$", x$`Ramsar Site No.`)), !anyDuplicated(x$`Ramsar Site No.`))
ok_date <- function(v) all(grepl("^\\d{4}-\\d{2}-\\d{2}$", v))
stopifnot(ok_date(x$`Designation date`), ok_date(x$`Last publication date`))
num <- function(v) { z <- suppressWarnings(as.numeric(v)); stopifnot(!anyNA(z[nzchar(v)])); z }
yes <- function(v) { stopifnot(all(v %in% c("Yes", "No"))); v == "Yes" }

# --- site table ----------------------------------------------------------------
s <- data.table(
  site_id = as.integer(x$`Ramsar Site No.`),
  site_name = trimws(x$`Site name`),
  region = x$Region, country = x$Country, territory = x$Territory,
  designation_date = as.IDate(x$`Designation date`),
  last_publication_date = as.IDate(x$`Last publication date`),
  area_ha = num(x$`Area (ha)`), lat = num(x$Latitude), lon = num(x$Longitude),
  max_elevation_m = num(x$`Maximum elevation`), min_elevation_m = num(x$`Minimum elevation`))
for (k in 1:9) set(s, j = paste0("criterion_", k), value = yes(x[[paste0("Criterion", k)]]))
s[, n_criteria := rowSums(.SD), .SDcols = paste0("criterion_", 1:9)]
# A few sites (mostly with sheets last published in the 1990s) have "No" for all
# nine criteria: the criteria were not recorded, which is not the same as none met.
s[, criteria_recorded := n_criteria > 0]
s[, montreux := yes(x$`Montreux listed`)]
s[, management_plan_implemented := yes(x$`Management plan implemented`)]
stopifnot(all(x$`Management plan available` %in% c("Yes", "No", "In preparation")))
s[, management_plan_available := x$`Management plan available`]
s[, transboundary := yes(x$`Does the wetland extend onto the territory of one or more other countries?`)]
s[, ramsar_advisory_mission := yes(x$`Ramsar Advisory Mission?`)]

# Other designations: "Overlapping site categories, <category heading>, item, item, ..."
items <- function(v, heading) lapply(strsplit(sub("^Overlapping site categories, ", "", v), ", ", fixed = TRUE),
                                     function(z) setdiff(z, heading))
glob <- items(x$`Global international legal designations`, "Global international legal and other formal")
regn <- items(x$`Regional international legal designations`, "Regional international legal and other formal designations")
natl <- items(x$`National conservation designation`, "Non-statutory designations")
has <- function(l, what) vapply(l, function(z) what %in% z, logical(1))
s[, `:=`(world_heritage = has(glob, "World Heritage site"), biosphere_reserve = has(glob, "UNESCO Biosphere Reserve"),
         other_global_designation = has(glob, "Other global designation"),
         natura2000 = has(regn, "EU Natura 2000"), other_international_designation = has(regn, "Other international designation"),
         important_bird_area = has(natl, "Important Bird Area"), important_plant_area = has(natl, "Important Plant Area"),
         other_nonstatutory_designation = has(natl, "Other non-statutory designation"))]
known <- c("World Heritage site", "UNESCO Biosphere Reserve", "Other global designation", "EU Natura 2000",
           "Other international designation", "Important Bird Area", "Important Plant Area", "Other non-statutory designation")
unknown <- setdiff(unique(c(unlist(glob), unlist(regn), unlist(natl))), c(known, ""))
if (length(unknown)) stop("Unknown designation items: ", paste(unknown, collapse = "; "))

# --- wetland types ---------------------------------------------------------------
# "Marine or coastal wetlands, A: Permanent shallow marine waters, G: ..., Inland wetlands, Ss: ..."
groups <- c("Marine or coastal wetlands" = "marine", "Inland wetlands" = "inland", "Human-made wetlands" = "human_made")
code_re <- "(?:[A-Z][a-z]?|[1-9]|Zk\\([abc]\\))"
valid <- list(marine = c(LETTERS[1:11], "Ga", "Zk(a)"),
              inland = c("L", "M", "N", "O", "P", "Q", "R", "Sp", "Ss", "Tp", "Ts", "U", "Va", "Vt", "W", "Xf", "Xp", "Y", "Zg", "Zk(b)"),
              human_made = c(as.character(1:9), "Zk(c)"))
split_re <- paste0(", (?=", code_re, ": |(?:Marine or coastal|Inland|Human-made) wetlands(?:, |$))")
wt <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
  v <- x$`Wetland Type`[i]
  if (!nzchar(v)) return(NULL)
  tok <- strsplit(v, split_re, perl = TRUE)[[1]]
  grp <- NA_character_; out <- list()
  for (t in tok) {
    if (t %in% names(groups)) { grp <- groups[[t]]; next }
    code <- sub(":.*", "", t); label <- sub("^[^:]+: ", "", t)
    if (is.na(grp) || !code %in% valid[[grp]]) stop(sprintf("Site %s: unexpected wetland type '%s'", x$`Ramsar Site No.`[i], t))
    out[[length(out) + 1]] <- data.table(site_id = as.integer(x$`Ramsar Site No.`[i]), group = grp, code = code, label = label)
  }
  rbindlist(out)
}))
wt <- unique(wt)
w <- wt[, .(n_wetland_types = .N, wetland_codes = paste(code, collapse = ";"),
            marine = any(group == "marine"), inland = any(group == "inland"), human_made = any(group == "human_made")),
        by = site_id]
s <- merge(s, w, by = "site_id", all.x = TRUE)
s[is.na(n_wetland_types), `:=`(n_wetland_types = 0L, marine = FALSE, inland = FALSE, human_made = FALSE)]
setorder(s, site_id)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fwrite(s, file.path(out_dir, "rsis_sites.csv"))
fwrite(wt[order(site_id)], file.path(out_dir, "rsis_wetland_types.csv"))

cat(sprintf("Sites: %d (numbers %d-%d), %d countries; designated %s to %s; sheets last published %s to %s\n",
            nrow(s), min(s$site_id), max(s$site_id), uniqueN(s$country), min(s$designation_date), max(s$designation_date),
            min(s$last_publication_date), max(s$last_publication_date)))
cat(sprintf("Criteria (sites): %s\n", paste(sprintf("%d: %d", 1:9, colSums(s[, paste0("criterion_", 1:9), with = FALSE])), collapse = ", ")))
cat(sprintf("No criterion recorded: %d; no wetland type: %d; Montreux Record: %d; management plan implemented: %d\n",
            sum(s$n_criteria == 0), sum(s$n_wetland_types == 0), sum(s$montreux), sum(s$management_plan_implemented)))
cat(sprintf("Wetland types: %d site x type rows; sites with marine %d, inland %d, human-made %d\n",
            nrow(wt), sum(s$marine), sum(s$inland), sum(s$human_made)))
