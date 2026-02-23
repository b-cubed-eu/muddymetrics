library(shiny)
library(leaflet)
library(dplyr)
library(ggplot2)
library(viridis)
library(plotly)
library(DT)
library(bslib)
library(sf)

ui <- fluidPage(
  theme = bslib::bs_theme(version = 5, bootswatch = "cosmo"),
  
  titlePanel("Muddymetrics: Ramsar Site Biodiversity Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Filters"),
      
      selectInput("continent",
                  "Select Continent:",
                  choices = c("All", "Africa", "Antarctica", "Asia", "Europe",
                              "North America", "Oceania", "South America"),
                  selected = "All"),
      
      selectInput("country",
                  "Select Country:",
                  choices = NULL),
      
      selectInput("site",
                  "Select Site:",
                  choices = NULL),
      
      hr(),
      
      h4("Summary"),
      verbatimTextOutput("summary_stats"),
      
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Global Map",
                 leafletOutput("global_map", height = 600)),
        
        tabPanel("Data Sufficiency",
                 plotlyOutput("density_hist", height = 400),
                 DT::dataTableOutput("sufficiency_table")),
        
        tabPanel("Site Details",
                 uiOutput("site_header"),
                 plotlyOutput("richness_ts", height = 300),
                 plotlyOutput("occurrence_ts", height = 300)),
        
        tabPanel("Country Comparison",
                 plotlyOutput("country_violin", height = 500))
      )
    )
  )
)


server <- function(input, output, session) {

  get_base_dir <- function() {
    app_dir <- getwd()
    if (basename(app_dir) == "dashboard") {
      return(file.path(app_dir, ".."))
    }
    return(app_dir)
  }

  get_site_data <- reactive({
    base_dir <- file.path(get_base_dir(), "inst", "extdata")
    shapefile_dir <- file.path(base_dir, "ramsar_sites_wkt")
    
    continent_map <- c(
      "ramsar_site_data_100m_africa" = "Africa",
      "ramsar_site_data_100m_antarctica" = "Antarctica",
      "ramsar_site_data_100m_asia" = "Asia",
      "ramsar_site_data_100m_europe" = "Europe",
      "ramsar_site_data_100m_northamerica" = "North America",
      "ramsar_site_data_100m_oceania" = "Oceania",
      "ramsar_site_data_100m_southamerica" = "South America"
    )
    
    sites <- list()
    
    for (data_dir_name in names(continent_map)) {
      data_dir <- file.path(base_dir, data_dir_name)
      continent <- continent_map[[data_dir_name]]
      
      if (!dir.exists(data_dir)) next
      
      countries <- list.files(data_dir)
      
      for (country in countries) {
        country_dir <- file.path(data_dir, country)
        if (!dir.exists(country_dir)) next
        
        files <- list.files(country_dir, pattern = "\\.csv$")
        
        for (f in files) {
          site_name <- sub("_data\\.csv$", "", f)
          site_id <- paste0(country, "_", site_name)
          
          wkt_file <- file.path(shapefile_dir, country, paste0(site_name, ".wkt"))
          
          coords <- NULL
          if (file.exists(wkt_file)) {
            tryCatch({
              wkt <- readLines(wkt_file, warn = FALSE)
              geom <- sf::st_as_sfc(wkt)
              centroids <- sf::st_centroid(geom)
              coords <- sf::st_coordinates(centroids)
            }, error = function(e) NULL)
          }
          
          sites[[length(sites) + 1]] <- list(
            site_id = site_id,
            site_name = site_name,
            country = country,
            continent = continent,
            lon = if (!is.null(coords)) coords[1] else NA,
            lat = if (!is.null(coords)) coords[2] else NA
          )
        }
      }
    }
    
    if (length(sites) == 0) {
      return(data.frame(
        site_id = character(),
        site_name = character(),
        country = character(),
        continent = character(),
        lon = numeric(),
        lat = numeric(),
        stringsAsFactors = FALSE
      ))
    }
    
    bind_rows(sites)
  })
  
  output$summary_stats <- renderText({
    data <- get_site_data()
    req(nrow(data) > 0)
    
    paste0(
      "Total Sites: ", nrow(data), "\n",
      "Countries: ", length(unique(data$country))
    )
  })
  
  observe({
    data <- get_site_data()
    continents <- unique(data$continent)
    updateSelectInput(session, "continent",
                     choices = c("All", sort(continents)))
  })
  
  observe({
    data <- get_site_data()
    cont <- input$continent
    
    if (cont == "All") {
      countries <- unique(data$country)
    } else {
      countries <- unique(data$country[data$continent == cont])
    }
    
    updateSelectInput(session, "country",
                     choices = c("All", sort(countries)))
  })
  
  observe({
    data <- get_site_data()
    cont <- input$continent
    country <- input$country
    
    if (cont == "All" && country == "All") {
      sites <- unique(data$site_id)
    } else if (cont != "All" && country == "All") {
      sites <- unique(data$site_id[data$continent == cont])
    } else if (cont == "All" && country != "All") {
      sites <- unique(data$site_id[data$country == country])
    } else {
      sites <- unique(data$site_id[data$continent == cont & data$country == country])
    }
    
    updateSelectInput(session, "site", choices = sort(sites))
  })
  
  output$global_map <- renderLeaflet({
    data <- get_site_data()
    
    if (input$continent != "All") {
      data <- data[data$continent == input$continent, ]
    }
    
    if (input$country != "All") {
      data <- data[data$country == input$country, ]
    }
    
    data <- data[!is.na(data$lon) & !is.na(data$lat), ]
    
    leaflet(data) |>
      addTiles() |>
      addCircleMarkers(
        ~lon, ~lat,
        popup = ~paste0(
          "<b>", site_name, "</b><br>",
          "Country: ", country
        ),
        radius = 5,
        color = "blue",
        fillOpacity = 0.7
      ) |>
      setView(lng = 10, lat = 50, zoom = 4)
  })
  
  output$density_hist <- renderPlotly({
    data <- get_site_data()
    
    plot_ly(data, x = ~country, type = "histogram") |>
      layout(
        title = "Sites by Country",
        xaxis = list(title = "Country"),
        yaxis = list(title = "Count")
      )
  })
  
  output$sufficiency_table <- DT::renderDataTable({
    data <- get_site_data()
    
    if (input$continent != "All") {
      data <- data[data$continent == input$continent, ]
    }
    
    if (input$country != "All") {
      data <- data[data$country == input$country, ]
    }
    
    data |>
      select(site_name, country, continent) |>
      DT::datatable()
  })
  
  output$site_header <- renderUI({
    data <- get_site_data()
    site_data <- data[data$site_id == input$site, ]
    
    if (nrow(site_data) == 0) return(NULL)
    
    wellPanel(
      h3(site_data$site_name),
      p(strong("Country:"), site_data$country),
      p(strong("Continent:"), site_data$continent),
      p(strong("Coordinates:"), 
        if (!is.na(site_data$lon)) paste0(
          round(site_data$lat, 4), ", ", round(site_data$lon, 4)
        ) else "N/A")
    )
  })
  
  find_site_csv <- function(site_id) {
    base_dir <- file.path(get_base_dir(), "inst", "extdata")
    
    continent_map <- c(
      "ramsar_site_data_100m_africa" = "Africa",
      "ramsar_site_data_100m_antarctica" = "Antarctica",
      "ramsar_site_data_100m_asia" = "Asia",
      "ramsar_site_data_100m_europe" = "Europe",
      "ramsar_site_data_100m_northamerica" = "North America",
      "ramsar_site_data_100m_oceania" = "Oceania",
      "ramsar_site_data_100m_southamerica" = "South America"
    )
    
    parts <- strsplit(site_id, "_")[[1]]
    country <- parts[1]
    site_name <- paste(parts[-1], collapse = "_")
    
    for (dir_name in names(continent_map)) {
      csv_file <- file.path(base_dir, dir_name, country, paste0(site_name, "_data.csv"))
      if (file.exists(csv_file)) {
        return(csv_file)
      }
    }
    return(NULL)
  }
  
  output$richness_ts <- renderPlotly({
    req(input$site)
    
    csv_file <- find_site_csv(input$site)
    
    if (!is.null(csv_file)) {
      site_data <- read.csv(csv_file)
      
      yearly <- site_data |>
        group_by(year) |>
        summarise(n_species = n_distinct(specieskey))
      
      plot_ly(yearly, x = ~year, y = ~n_species, type = "scatter", mode = "lines+markers") |>
        layout(
          title = "Species Richness Over Time",
          xaxis = list(title = "Year"),
          yaxis = list(title = "Number of Species")
        )
    } else {
      plot_ly() |>
        add_text(x = 0.5, y = 0.5, text = "No data available") |>
        layout(title = "Species Richness Over Time")
    }
  })
  
  output$occurrence_ts <- renderPlotly({
    req(input$site)
    
    csv_file <- find_site_csv(input$site)
    
    if (!is.null(csv_file)) {
      site_data <- read.csv(csv_file)
      
      yearly <- site_data |>
        group_by(year) |>
        summarise(total_occ = sum(occurrences))
      
      plot_ly(yearly, x = ~year, y = ~total_occ, type = "scatter", mode = "lines+markers") |>
        layout(
          title = "Total Occurrences Over Time",
          xaxis = list(title = "Year"),
          yaxis = list(title = "Total Occurrences")
        )
    } else {
      plot_ly() |>
        add_text(x = 0.5, y = 0.5, text = "No data available") |>
        layout(title = "Total Occurrences Over Time")
    }
  })
  
  output$country_violin <- renderPlotly({
    data <- get_site_data()
    
    country_counts <- data |>
      group_by(country) |>
      summarise(n = n())
    
    plot_ly(country_counts, x = ~reorder(country, -n), y = ~n, type = "bar") |>
      layout(
        title = "Number of Sites by Country",
        xaxis = list(title = "Country", tickangle = 45),
        yaxis = list(title = "Number of Sites")
      )
  })
}


shinyApp(ui = ui, server = server)
