#!/usr/bin/env Rscript

# Generate static HTML gallery from pre-generated PNG plots

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

get_site_name <- function(filename) {
  sub("_data.*$", "", sub(".*site_[0-9]+_", "", filename))
}

get_plot_type <- function(filename) {
  if (grepl("obs_richness", filename)) return("Observed Richness")
  if (grepl("cum_richness", filename)) return("Cumulative Richness")
  if (grepl("total_occ_ts", filename)) return("Total Occurrences (Time Series)")
  if (grepl("total_occ_map", filename)) return("Total Occurrences (Map)")
  if (grepl("overall_density", filename)) return("Overall Density")
  return("Other")
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

cat("Generating index page...\n")

index_content <- html_header("Muddymetrics - Ramsar Site Biodiversity Dashboard")

index_content <- paste0(index_content, '
  <header>
    <h1>🧪 Muddymetrics</h1>
    <p>Ramsar Site Biodiversity Indicators</p>
  </header>
  <nav>
    <a href="index.html">Home</a>
  </nav>
  <main>
    <div class="stats">
      <div class="stats-grid>
')

total_sites <- 0
for (continent in continents) {
  continent_dir <- file.path(output_base, continent)
  if (!dir.exists(continent_dir)) next
  
  countries <- list.files(continent_dir)
  for (country in countries) {
    country_dir <- file.path(continent_dir, country)
    if (!dir.exists(country_dir)) next
    files <- list.files(country_dir, pattern = "\\.png$")
    sites <- unique(gsub("_.*$", "", gsub("site_[0-9]+_", "", files)))
    total_sites <- total_sites + length(sites)
  }
}

index_content <- paste0(index_content, '
        <div class="stat">
          <div class="number">', total_sites, '</div>
          <div class="label">Sites</div>
        </div>
        <div class="stat">
          <div class="number">7</div>
          <div class="label">Continents</div>
        </div>
        <div class="stat">
          <div class="number">~2000+</div>
          <div class="label">Plots</div>
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

cat("Generating continent pages...\n")

for (continent in continents) {
  continent_dir <- file.path(output_base, continent)
  if (!dir.exists(continent_dir)) next
  
  continent_html_dir <- file.path(docs_dir, continent)
  if (!dir.exists(continent_html_dir)) {
    dir.create(continent_html_dir, recursive = TRUE)
  }
  
  countries <- sort(list.files(continent_dir))
  
  continent_content <- html_header(paste0(continent_names[continent], " - Muddymetrics"))
  continent_content <- paste0(continent_content, '
  <header>
    <h1>', continent_names[continent], '</h1>
  </header>
  <nav>
    <a href="../index.html">← Home</a>
  </nav>
  <main>
    <div class="breadcrumb">
      <a href="../index.html">Home</a> / ', continent_names[continent], '
    </div>
    <h2>Countries in ', continent_names[continent], '</h2>
    <div class="grid">
  ')
  
  for (country in countries) {
    country_dir <- file.path(continent_dir, country)
    if (!dir.exists(country_dir)) next
    
    png_files <- list.files(country_dir, pattern = "\\.png$")
    n_plots <- length(png_files)
    
    if (n_plots == 0) next
    
    continent_content <- paste0(continent_content, '
      <div class="card">
        <a href="', country, '/index.html">
          <div class="card-content">
            <h3>', country, '</h3>
            <span class="type">', n_plots, ' plots</span>
          </div>
        </a>
      </div>
    ')
  }
  
  continent_content <- paste0(continent_content, '
    </div>
  </main>
  ', html_footer())
  
  writeLines(continent_content, file.path(continent_html_dir, "index.html"))
  
  cat("  ", continent, ":", length(countries), "countries\n")
  
  for (country in countries) {
    country_dir <- file.path(continent_dir, country)
    if (!dir.exists(country_dir)) next
    
    country_html_dir <- file.path(continent_html_dir, country)
    if (!dir.exists(country_html_dir)) {
      dir.create(country_html_dir, recursive = TRUE)
    }
    
    png_files <- list.files(country_dir, pattern = "\\.png$")
    
    if (length(png_files) == 0) next
    
    country_content <- html_header(paste0(country, " - ", continent_names[continent], " - Muddymetrics"))
    country_content <- paste0(country_content, '
  <header>
    <h1>', country, '</h1>
    <p>', continent_names[continent], '</p>
  </header>
  <nav>
    <a href="../index.html">← ', continent_names[continent], '</a>
    <a href="../../index.html">← Home</a>
  </nav>
  <main>
    <div class="breadcrumb">
      <a href="../../index.html">Home</a> / 
      <a href="../index.html">', continent_names[continent], '</a> / 
      ', country, '
    </div>
    <h2>Site Plots in ', country, '</h2>
    <div class="grid">
  ')
    
    site_ids <- unique(gsub("_data.*$", "", sub("site_([0-9]+)_.*", "site_\\1", png_files)))
    
    for (site_id in site_ids) {
      site_pngs <- grep(paste0("^", site_id, "_"), png_files, value = TRUE)
      
      if (length(site_pngs) == 0) next
      
      site_name <- gsub("_", " ", sub("site_[0-9]+_", "", site_id))
      
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
    
    for (site_id in site_ids) {
      site_pngs <- grep(paste0("^", site_id, "_"), png_files, value = TRUE)
      
      if (length(site_pngs) == 0) next
      
      site_html_dir <- file.path(country_html_dir, site_id)
      if (!dir.exists(site_html_dir)) {
        dir.create(site_html_dir, recursive = TRUE)
      }
      
      site_name <- gsub("_", " ", sub("site_[0-9]+_", "", site_id))
      
      site_content <- html_header(paste0(site_name, " - ", country, " - Muddymetrics"))
      site_content <- paste0(site_content, '
  <header>
    <h1>', site_name, '</h1>
    <p>', country, ', ', continent_names[continent], '</p>
  </header>
  <nav>
    <a href="../index.html">← ', country, '</a>
    <a href="../../index.html">← Home</a>
  </nav>
  <main>
    <div class="breadcrumb">
      <a href="../../index.html">Home</a> / 
      <a href="../index.html">', continent_names[continent], '</a> / 
      <a href="../', country, '/index.html">', country, '</a> / 
      ', site_name, '
    </div>
    <h2>Plots for ', site_name, '</h2>
    <div class="grid">
  ')
      
      for (png in sort(site_pngs)) {
        plot_type <- get_plot_type(png)
        
        site_content <- paste0(site_content, '
      <div class="card">
        <a href="../', png, '" target="_blank">
          <div class="card-content">
            <h3>', plot_type, '</h3>
            <span class="type">View Plot</span>
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
}

cat("\n✓ Gallery generated in docs/\n")
cat("To preview locally: serve the docs/ folder\n")
cat("To deployHub and enable Pages\n")
: push to Git