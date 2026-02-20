library(shiny)
library(leaflet)
library(dplyr)
library(ggplot2)
library(viridis)
library(plotly)

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

      h4("Data Sufficiency"),
      verbatimTextOutput("site_summary"),

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
                 plotlyOutput("country_violin", height = 500)),

        tabPanel("About",
                 includeMarkdown("README.md"))
      )
    )
  )
)


server <- function(input, output, session) {

  summary_data <- reactiveFileReader(
    intervalMillis = 5000,
    session = session,
    filePath = "output/global_sufficiency_summary.csv",
    readFunc = read.csv
  )

  observe({
    data <- summary_data()
    continents <- unique(data$continent)
    updateSelectInput(session, "continent",
                     choices = c("All", sort(continents)))
  })

  observe({
    data <- summary_data()
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
    data <- summary_data()
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
    data <- summary_data()

    if (input$continent != "All") {
      data <- data[data$continent == input$continent, ]
    }

    if (input$country != "All") {
      data <- data[data$country == input$country, ]
    }

    pal <- colorFactor(c("red", "green"), domain = c("Data-Poor", "Data-Rich"))

    leaflet(data) |>
      addTiles() |>
      addCircleMarkers(
        ~as.numeric(sub(".*_", "", sub("site_", "", site_id))),
        ~as.numeric(sub(".*_", "", sub("site_", "", site_id))) * 0.5,
        color = ~pal(data_class),
        popup = ~paste0(
          "<b>", site_name, "</b><br>",
          "Country: ", country, "<br>",
          "Density: ", round(density_km2, 3), "<br>",
          "Status: ", data_class
        ),
        radius = 5
      ) |>
      addLegend(position = "bottomright",
                pal = pal,
                values = ~data_class,
                title = "Data Status")
  })

  output$density_hist <- renderPlotly({
    data <- summary_data()

    if (input$continent != "All") {
      data <- data[data$continent == input$continent, ]
    }

    if (input$country != "All") {
      data <- data[data$country == input$country, ]
    }

    plot_ly(data, x = ~density_km2, type = "histogram",
            marker = list(color = "#1f77b4")) |>
      layout(
        title = "Distribution of Occurrence Density",
        xaxis = list(title = "Density (records/km²)"),
        yaxis = list(title = "Count")
      )
  })

  output$sufficiency_table <- DT::renderDataTable({
    data <- summary_data()

    if (input$continent != "All") {
      data <- data[data$continent == input$continent, ]
    }

    if (input$country != "All") {
      data <- data[data$country == input$country, ]
    }

    data |>
      select(site_name, country, density_km2, data_class) |>
      DT::datatable()
  })

  output$site_header <- renderUI({
    data <- summary_data()
    site_data <- data[data$site_id == input$site, ]

    if (nrow(site_data) == 0) return(NULL)

    wellPanel(
      h3(site_data$site_name),
      p(strong("Country:"), site_data$country),
      p(strong("Density:"), round(site_data$density_km2, 3), "records/km²"),
      p(strong("Status:"), site_data$data_class)
    )
  })

  output$richness_ts <- renderPlotly({
    plot_ly(x = 1:10, y = rnorm(10), type = "scatter", mode = "lines") |>
      layout(title = "Species Richness Over Time",
             xaxis = list(title = "Year"),
             yaxis = list(title = "Species Richness"))
  })

  output$occurrence_ts <- renderPlotly({
    plot_ly(x = 1:10, y = rnorm(10), type = "scatter", mode = "lines") |>
      layout(title = "Total Occurrences Over Time",
             xaxis = list(title = "Year"),
             yaxis = list(title = "Occurrences"))
  })

  output$country_violin <- renderPlotly({
    data <- summary_data()

    if (input$continent != "All") {
      data <- data[data$continent == input$continent, ]
    }

    plot_ly(data, x = ~country, y = ~density_km2, type = "violin") |>
      layout(title = "Density by Country",
             xaxis = list(title = "Country"),
             yaxis = list(title = "Density (records/km²)"))
  })
}


shinyApp(ui = ui, server = server)
