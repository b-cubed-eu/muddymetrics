# GDP per capita (current US$, World Bank indicator NY.GDP.PCAP.CD) for the
# countries of the Ramsar Sites, for the model of which sites can be verified.
#
# Usage (from the repository folder):
#   Rscript scripts/get_worldbank_gdp.R [--countries _local_data/RSIS/rsis_country_iso2.csv]
#           [--out _local_data/gbif_country/worldbank_gdp_per_capita.csv] [--years 2019:2023]
# Needs jsonlite. Keeps the raw API response next to the output. For each
# country the latest year with a value within --years is used (gdp_year says which).

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Please install the jsonlite package")
args <- commandArgs(TRUE)
opt <- function(name, default) { i <- match(paste0("--", name), args); if (is.na(i)) default else args[i + 1] }
COUNTRIES <- opt("countries", "_local_data/RSIS/rsis_country_iso2.csv")
OUT <- opt("out", "_local_data/gbif_country/worldbank_gdp_per_capita.csv")
YEARS <- opt("years", "2019:2023")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)

raw <- sub("\\.csv$", "_raw.json", OUT)
url <- sprintf("https://api.worldbank.org/v2/country/all/indicator/NY.GDP.PCAP.CD?format=json&date=%s&per_page=20000", YEARS)
if (!file.exists(raw)) utils::download.file(url, raw, quiet = TRUE, mode = "wb", method = "libcurl")
j <- jsonlite::fromJSON(raw)
stopifnot(length(j) == 2, j[[1]]$pages == 1)
d <- j[[2]]
wb <- data.frame(iso2 = d$country$id, wb_name = d$country$value, year = as.integer(d$date), gdp = d$value, stringsAsFactors = FALSE)
wb <- wb[!is.na(wb$gdp), ]
wb <- wb[order(wb$iso2, -wb$year), ]
wb <- wb[!duplicated(wb$iso2), ]

# na.strings = "": Namibia's code is "NA".
cn <- utils::read.csv(COUNTRIES, encoding = "UTF-8", stringsAsFactors = FALSE, na.strings = "")
out <- merge(cn, wb, by = "iso2", all.x = TRUE)
names(out)[names(out) == "gdp"] <- "gdp_per_capita_usd"; names(out)[names(out) == "year"] <- "gdp_year"
out$indicator <- "NY.GDP.PCAP.CD"; out$retrieved <- format(Sys.Date())
utils::write.csv(out[order(out$country), ], OUT, row.names = FALSE, fileEncoding = "UTF-8")
cat(sprintf("Countries: %d; with GDP per capita: %d; latest year used: %s\n", nrow(out), sum(!is.na(out$gdp_per_capita_usd)),
            paste(names(table(out$gdp_year)), table(out$gdp_year), sep = ": ", collapse = ", ")))
if (any(is.na(out$gdp_per_capita_usd))) cat("Without a value:", paste(out$country[is.na(out$gdp_per_capita_usd)], collapse = "; "), "\n")
