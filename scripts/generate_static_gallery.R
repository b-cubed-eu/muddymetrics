#!/usr/bin/env Rscript

# Generate static HTML gallery from PNG plots

output_base <- "output/ramsar_metric_results_100m"
docs_dir <- "docs"

if (!dir.exists(docs_dir)) {
  dir.create(docs_dir)
}

continents <- c("africa", "antarctica", "asia", "europe", "northamerica", "oceania", "southamerica")

continent_names <- setNames(
  c("Africa", "Antarctica", "Asia", "Europe", "North America", "Oceania", "South America"),
  continents
)

get_plot_type <- function(filename) {
  if (grepl("obs_richness", filename)) return("Observed Richness")
  if (grepl("cum_richness", filename)) return("Cumulative Richness")
  if (grepl("total_occ_ts", filename)) return("Total Occurrences")
  if (grepl("total_occ_map", filename)) return("Total Occurrences Map")
  if (grepl("overall_density", filename)) return("Overall Density")
  if (grepl("completeness", filename)) return("Completeness")
  return("Other")
}

extract_site_id <- function(filename) {
  # Extract site_XXXX from filename like site_1001_Cromarty_Firth_total_occ_ts.png
  if (grepl("^site_[0-9]+_", filename)) {
    return(sub("^(site_[0-9]+)_.*", "\\1", filename))
  }
  return(NULL)
}

html_header <- function(title) {
  paste0('<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>', title, '</title>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; 
           background: #f5f5f5; color: #333; line-height: 1.6; }
    header { background: #1a1a2e; color: white; padding: 2rem; text-align: center; }
    header h1 { font-size: 2rem; margin-bottom: 0.5rem; }
    header p { color: #aaa; }
    nav { background: #16213e; padding: 1rem; text-align: center; }
    nav a { color: white; text-decoration: none; margin: 0 1rem; padding: 0.5rem 1rem; 
            background: #0f3460; border-radius: 4px; }
    nav a:hover { background: #e94560; }
    main { max-width: 1200px; margin: 2rem auto; padding: 0 1rem; }
    h2 { margin-bottom: 1rem; color: #1a1a2e; border-bottom: 2px solid #e94560; padding-bottom: 0.5rem; }
    .breadcrumb { padding: 1rem; background: white; border-radius: 8px; margin-bottom: 1rem; }
    .breadcrumb a { color: #0f3460; text-decoration: none; }
    .breadcrumb a:hover { text-decoration: underline; }
    .grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(300px, 1fr)); gap: 1.5rem; }
    .card { background: white; border-radius: 8px; overflow: hidden; box-shadow: 0 2px 8px rgba(0,0,0,0.1); 
            transition: transform 0.2s; }
    .card:hover { transform: translateY(-4px); box-shadow: 0 4px 16px rgba(0,0,0,0.15); }
    .card img { width: 100%; height: 200px; object-fit: cover; }
    .card-content { padding: 1rem; }
    .card h3 { font-size: 1rem; margin-bottom: 0.5rem; color: #1a1a2e; }
    .card .type { display: inline-block; background: #e94560; color: white; 
                  padding: 0.25rem 0.5rem; border-radius: 4px; font-size: 0.75rem; }
    .card a { text-decoration: none; color: inherit; }
    footer { text-align: center; padding: 2rem; color: #666; margin-top: 2rem; }
    .stats { background: white; padding: 1rem; border-radius: 8px; margin-bottom: 2rem; }
    .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 1rem; }
    .stat { text-align: center; }
    .stat .number { font-size: 2rem; font-weight: bold; color: #e94560; }
    .stat .label { color: #666; font-size: 0.875rem; }
  </style>
</head>
<body>
')
}

html_footer <- function() {
  '
  <footer>
    <p>Muddymetrics - Ramsar Site Biodiversity Dashboard</p>
  </footer>
</body>
</html>
'
}

cat("Generating static gallery...\n")

# Generate index
total_sites <- 0
for (continent in continents) {
  continent_dir <- file.path(output_base, continent)
  if (!dir.exists(continent_dir)) next
  countries <- list.files(continent_dir)
  for (country in countries) {
    country_dir <- file.path(continent_dir, country)
    if (!dir.exists(country_dir)) next
    png_files <- list.files(country_dir, pattern = "\\.png$")
    site_ids <- unique(sapply(png_files, extract_site_id))
    site_ids <- site_ids[!is.null(site_ids)]
    total_sites <- total_sites + length(site_ids)
  }
}

cat("Total sites:", total_sites, "\n")

# Generate index page
index_content <- html_header("Muddymetrics - Ramsar Site Biodiversity Dashboard")
index_content <- paste0(index_content, '
  <header>
    <h1>Muddymetrics</h1>
    <p>Ramsar Site Biodiversity Indicators</p>
  </header>
  <nav>
    <a href="index.html">Home</a>
  </nav>
  <main>
    <div class="stats">
      <div class="stats-grid">
        <div class="stat">
          <div class="number">', total_sites, '</div>
          <div class="label">Sites</div>
        </div>
        <div class="stat">
          <div class="number">7</div>
          <div class="label">Continents</div>
        </div>
      </div>
    </div>
    <h2>Browse by Continent</h2>
    <div class="grid">
')

for (continent in continents) {
  continent_dir <- file.path(output_base, continent)
  if (!dir.exists(continent_dir)) next
  countries <- list.files(continent_dir)
  n_countries <- length(countries)
  index_content <- paste0(index_content, '
      <div class="card">
        <a href="', continent, '/index.html">
          <div class="card-content">
            <h3>', continent_names[continent], '</h3>
            <span class="type">', n_countries, ' countries</span>
          </div>
        </a>
      </div>
  ')
}

index_content <- paste0(index_content, '
    </div>
  </main>
', html_footer())
writeLines(index_content, file.path(docs_dir, "index.html"))

# Process continents
for (continent in continents) {
  continent_dir <- file.path(output_base, continent)
  if (!dir.exists(continent_dir)) next
  
  continent_html_dir <- file.path(docs_dir, continent)
  if (!dir.exists(continent_html_dir)) {
    dir.create(continent_html_dir, recursive = TRUE)
  }
  
  countries <- sort(list.files(continent_dir))
  
  # Continent index
  cont_content <- html_header(paste0(continent_names[continent], " - Muddymetrics"))
  cont_content <- paste0(cont_content, '
  <header>
    <h1>', continent_names[continent], '</h1>
  </header>
  <nav>
    <a href="../index.html">Home</a>
  </nav>
  <main>
    <div class="breadcrumb">
      <a href="../index.html">Home</a> / ', continent_names[continent], '
    </div>
    <h2>Countries</h2>
    <div class="grid">
')
  
  for (country in countries) {
    country_dir <- file.path(continent_dir, country)
    if (!dir.exists(country_dir)) next
    png_files <- list.files(country_dir, pattern = "\\.png$")
    if (length(png_files) == 0) next
    
    # URL encode country name for links (not folder names)
    country_link <- URLencode(country)
    cont_content <- paste0(cont_content, '
      <div class="card">
        <a href="', country, '/index.html">
          <div class="card-content">
            <h3>', country, '</h3>
            <span class="type">', length(png_files), ' plots</span>
          </div>
        </a>
      </div>
    ')
  }
  
  cont_content <- paste0(cont_content, '
    </div>
  </main>
  ', html_footer())
  writeLines(cont_content, file.path(continent_html_dir, "index.html"))
  
  # Countries
  for (country in countries) {
    country_dir <- file.path(continent_dir, country)
    if (!dir.exists(country_dir)) next
    png_files <- list.files(country_dir, pattern = "\\.png$")
    if (length(png_files) == 0) next
    
    country_html_dir <- file.path(continent_html_dir, country)
    if (!dir.exists(country_html_dir)) {
      dir.create(country_html_dir, recursive = TRUE)
    }
    
    # Get unique sites
    site_ids <- unique(sapply(png_files, extract_site_id))
    site_ids <- site_ids[!is.null(site_ids)]
    
    # Country index
    country_content <- html_header(paste0(country, " - Muddymetrics"))
    country_content <- paste0(country_content, '
  <header>
    <h1>', country, '</h1>
    <p>', continent_names[continent], '</p>
  </header>
  <nav>
    <a href="../index.html">', continent_names[continent], '</a>
    <a href="../../index.html">Home</a>
  </nav>
  <main>
    <div class="breadcrumb">
      <a href="../../index.html">Home</a> / 
      <a href="../index.html">', continent_names[continent], '</a> / 
      ', country, '
    </div>
    <h2>Sites</h2>
    <div class="grid">
')
    
    for (site_id in site_ids) {
      site_pngs <- grep(paste0("^", site_id, "_"), png_files, value = TRUE)
      if (length(site_pngs) == 0) next
      
      site_name <- sub(paste0(site_id, "_"), "", site_id)
      site_name <- gsub("_", " ", site_name)
      
      country_content <- paste0(country_content, '
      <div class="card">
        <a href="', site_id, '/index.html">
          <div class="card-content">
            <h3>', site_name, '</h3>
            <span class="type">', length(site_pngs), ' plots</span>
          </div>
        </a>
      </div>
')
    }
    
    country_content <- paste0(country_content, '
    </div>
  </main>
  ', html_footer())
    writeLines(country_content, file.path(country_html_dir, "index.html"))
    
    # Site pages
    for (site_id in site_ids) {
      site_pngs <- grep(paste0("^", site_id, "_"), png_files, value = TRUE)
      if (length(site_pngs) == 0) next
      
      site_html_dir <- file.path(country_html_dir, site_id)
      if (!dir.exists(site_html_dir)) {
        dir.create(site_html_dir, recursive = TRUE)
      }
      
      # Copy PNG files
      for (png in site_pngs) {
        from_path <- file.path(country_dir, png)
        to_path <- file.path(site_html_dir, png)
        if (!file.exists(to_path)) {
          file.copy(from_path, to_path)
        }
      }
      
      site_name <- sub(paste0(site_id, "_"), "", site_id)
      site_name <- gsub("_", " ", site_name)
      
      site_content <- html_header(paste0(site_name, " - ", country))
      site_content <- paste0(site_content, '
  <header>
    <h1>', site_name, '</h1>
    <p>', country, ', ', continent_names[continent], '</p>
  </header>
  <nav>
    <a href="index.html">', country, '</a>
    <a href="../../index.html">Home</a>
  </nav>
  <main>
    <div class="breadcrumb">
      <a href="../../index.html">Home</a> / 
      <a href="index.html">', continent_names[continent], '</a> / 
      <a href="index.html">', country, '</a> / 
      ', site_name, '
    </div>
    <h2>Plots</h2>
    <div class="grid">
')
      
      for (png in sort(site_pngs)) {
        plot_type <- get_plot_type(png)
        site_content <- paste0(site_content, '
      <div class="card">
        <a href="', png, '" target="_blank">
          <img src="', png, '" alt="', plot_type, '">
          <div class="card-content">
            <h3>', plot_type, '</h3>
          </div>
        </a>
      </div>
')
      }
      
      site_content <- paste0(site_content, '
    </div>
  </main>
  ', html_footer())
      writeLines(site_content, file.path(site_html_dir, "index.html"))
    }
  }
  
  cat("  ", continent, ":", length(countries), "countries\n")
}

cat("\nDone! Serve with: cd docs && python -m http.server 8000\n")
