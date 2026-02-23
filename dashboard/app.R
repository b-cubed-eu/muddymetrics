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
      
      hr(),
      
      h4("Selected Site"),
      verbatimTextOutput("selected_site_info"),
      
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Map & Sites",
                 p("Click on a marker to select a site. Use the sidebar to filter by continent/country."),
                 leafletOutput("global_map", height = 500),
                 DT::dataTableOutput("sites_table", height = 300)),
        
        tabPanel("Site Details",
                 uiOutput("site_header"),
                 plotlyOutput("richness_ts", height = 300),
                 plotlyOutput("occurrence_ts", height = 300)),
        
        tabPanel("Country Overview",
                 plotlyOutput("sites_by_country", height = 400),
                 plotlyOutput("occurrences_by_country", height = 400))
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
  
  output$selected_site_info <- renderText({
    selected <- get_current_site_id()
    
    if (is.null(selected)) {
      return("No site selected\n\nClick on map\nor table to select")
    }
    
    data <- get_site_data()
    site_data <- data[data$site_id == selected, ]
    
    if (nrow(site_data) == 0) {
      return("No site selected\n\nClick on map\nor table to select")
    }
    
    paste0(
      "Site: ", site_data$site_name, "\n",
      "Country: ", site_data$country, "\n",
      "Continent: ", site_data$continent
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
  
  get_filtered_data <- reactive({
    data <- get_site_data()
    
    if (input$continent != "All") {
      data <- data[data$continent == input$continent, ]
    }
    
    if (input$country != "All") {
      data <- data[data$country == input$country, ]
    }
    
    data
  })
  
  output$global_map <- renderLeaflet({
    data <- get_filtered_data()
    data <- data[!is.na(data$lon) & !is.na(data$lat), ]
    
    leaflet(data) |>
      addTiles() |>
      addCircleMarkers(
        ~lon, ~lat,
        layerId = ~site_id,
        popup = ~paste0(
          "<b>", site_name, "</b><br>",
          "Country: ", country, "<br>",
          "<i>Click to select</i>"
        ),
        radius = 6,
        color = "#2E86AB",
        fillOpacity = 0.8
      ) |>
      setView(lng = 10, lat = 50, zoom = 4)
  })
  
  output$sites_table <- DT::renderDataTable({
    data <- get_filtered_data()
    
    data |>
      select(site_name, country, continent, lon, lat) |>
      mutate(lat = round(lat, 4), lon = round(lon, 4)) |>
      DT::datatable(
        rownames = FALSE,
        colnames = c("Site Name", "Country", "Continent", "Latitude", "Longitude"),
        selection = "single",
        server = FALSE,
        options = list(pageLength = 10)
      )
  })
  
  output$site_header <- renderUI({
    selected_site <- NULL
    
    if (!is.null(input$global_map_marker_click)) {
      selected_site <- input$global_map_marker_click$id
    } else if (!is.null(input$sites_table_rows_selected)) {
      data <- get_filtered_data()
      selected_site <- data$site_id[input$sites_table_rows_selected[1]]
    }
    
    if (is.null(selected_site)) {
      return(div(
        p("Select a site from the map or table above to view details."),
        style = "color: gray; padding: 20px;"
      ))
    }
    
    data <- get_site_data()
    site_data <- data[data$site_id == selected_site, ]
    
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
  
  get_current_site_id <- reactive({
    if (!is.null(input$global_map_marker_click)) {
      return(input$global_map_marker_click$id)
    }
    if (!is.null(input$sites_table_rows_selected)) {
      data <- get_filtered_data()
      return(data$site_id[input$sites_table_rows_selected[1]])
    }
    return(NULL)
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
    site_id <- get_current_site_id()
    req(site_id)
    
    csv_file <- find_site_csv(site_id)
    
    if (!is.null(csv_file)) {
      site_data <- read.csv(csv_file)
      
      yearly <- site_data |>
        group_by(year) |>
        summarise(n_species = n_distinct(specieskey))
      
      plot_ly(yearly, x = ~year, y = ~n_species, type = "scatter", mode = "lines+markers",
              marker = list(color = "#2E86AB"), line = list(color = "#2E86AB")) |>
        layout(
          title = "Species Richness Over Time",
          xaxis = list(title = "Year"),
          yaxis = list(title = "Number of Species"),
          margin = list(t = 40)
        )
    } else {
      plot_ly() |>
        add_text(x = 0.5, y = 0.5, text = "No data available", showlegend = FALSE) |>
        layout(title = "Species Richness Over Time")
    }
  })
  
  output$occurrence_ts <- renderPlotly({
    site_id <- get_current_site_id()
    req(site_id)
    
    csv_file <- find_site_csv(site_id)
    
    if (!is.null(csv_file)) {
      site_data <- read.csv(csv_file)
      
      yearly <- site_data |>
        group_by(year) |>
        summarise(total_occ = sum(occurrences))
      
      plot_ly(yearly, x = ~year, y = ~total_occ, type = "scatter", mode = "lines+markers",
              marker = list(color = "#A23B72"), line = list(color = "#A23B72")) |>
        layout(
          title = "Total Occurrences Over Time",
          xaxis = list(title = "Year"),
          yaxis = list(title = "Total Occurrences"),
          margin = list(t = 40)
        )
    } else {
      plot_ly() |>
        add_text(x = 0.5, y = 0.5, text = "No data available", showlegend = FALSE) |>
        layout(title = "Total Occurrences Over Time")
    }
  })
  
  output$sites_by_country <- renderPlotly({
    data <- get_filtered_data()
    
    country_counts <- data |>
      group_by(country) |>
      summarise(n = n()) |>
      arrange(desc(n))
    
    plot_ly(country_counts, x = ~reorder(country, -n), y = ~n, type = "bar",
            marker = list(color = "#2E86AB")) |>
      layout(
        title = "Number of Sites by Country",
        xaxis = list(title = "Country", tickangle = 45),
        yaxis = list(title = "Number of Sites"),
        margin = list(b = 100)
      )
  })
  
  output$occurrences_by_country <- renderPlotly({
    data <- get_filtered_data()
    
    country_occ <- list()
    
    for (i in 1:nrow(data)) {
      csv_file <- find_site_csv(data$site_id[i])
      if (!is.null(csv_file)) {
        tryCatch({
          site_data <- read.csv(csv_file)
          total_occ <- sum(site_data$occurrences, na.rm = TRUE)
          country_occ[[data$country[i]]] <- c(country_occ[[data$country[i]]], total_occ)
        }, error = function(e) NULL)
      }
    }
    
    occ_sum <- data.frame(
      country = names(country_occ),
      total_occ = sapply(country_occ, sum),
      stringsAsFactors = FALSE
    ) |>
      arrange(desc(total_occ))
    
    if (nrow(occ_sum) == 0) {
      return(plot_ly() |>
        add_text(x = 0.5, y = 0.5, text = "No occurrence data available") |>
        layout(title = "Total Occurrences by Country"))
    }
    
    plot_ly(occ_sum, x = ~reorder(country, -total_occ), y = ~total_occ, type = "bar",
            marker = list(color = "#A23B72")) |>
      layout(
        title = "Total Occurrences by Country",
        xaxis = list(title = "Country", tickangle = 45),
        yaxis = list(title = "Total Occurrences"),
        margin = list(b = 100)
      )
  })
}


shinyApp(ui = ui, server = server)
