# Load necessary libraries
# install.packages(c("tidyverse", "purrr", "ggplot2"), dependencies = TRUE)
library(tidyverse)
library(purrr)
library(ggplot2)

# --- 1. YOUR REAL DATA ASSIGNMENT ---
# Replace the code below with the assignment of your actual list.
#
# IMPORTANT: Your list MUST be structured as:
# occurrence_list <- list(
#   "Country A" = list("Site_1" = 10.5, "Site_2" = 8.1),
#   "Country B" = list("Site_3" = 22.0, "Site_4" = 19.5),
#   # ... and so on for all 30 countries and 200 sites
# )
#
# Please paste your actual list assignment here:
# ------------------------------------------------------------------------
# V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V V
set.seed(42)

# Create 30 mock country names
country_names <- paste("Country", sprintf("%02d", 1:30))

# Function to generate variable number of sites (between 3 and 10) for each country
generate_site_data <- function(country_name) {
  num_sites <- sample(3:10, 1)

  # Generate random occurrence values and append units to simulate your data format
  values <- map(1:num_sites, ~ {
    # 1. 20% chance of being NULL
    if (runif(1) < 0.2) return(NULL)

    # 2. Otherwise, return a value with units as a string
    val <- abs(rnorm(1, mean = runif(1, 5, 25), sd = runif(1, 1, 5)))
    paste0(val, " [1/km^2]")
  })

  # Name the sites Site_01, Site_02, etc.
  site_names <- paste0("site_", sprintf("%04d", 1:num_sites), "_data.csv")
  names(values) <- site_names

  # Return the inner list of sites
  return(values)
}

# Apply the function to create the final nested list
occurrence_list <- set_names(
  map(country_names, generate_site_data),
  country_names
)
# ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^
# ------------------------------------------------------------------------

create_violin_plot <- function(nested_list) {

  library(tidyverse)
  library(ggplot2)
  library(purrr)

  # Assign the input nested list to a variable for processing
  x <- nested_list

  # Check the structure of the input data
  cat("--- Input Data Structure ---\n")
  str(x)
  cat("\nTotal countries:", length(x), "\n\n")

  # --- 2. DATA TRANSFORMATION: Flatten the nested list and clean values ---

  # Get the list of country names from the input data
  country_names_from_list <- names(x)

  # Use map_dfr to iterate through the countries
  df_occurrences <- purrr::map_dfr(x, ~ {
    # .x is the inner list of sites for one country (e.g., $Vietnam)

    site_names <- names(.x)
    raw_values <- unname(.x)

    # 1. Process raw values: Handle NULL, strip units, and convert to numeric
    cleaned_values <- map_dbl(raw_values, function(val) {
      if (is.null(val) || is.na(val)) {
        return(NA_real_)
      }

      # Convert to string and use regular expression to remove the unit part
      # We specifically look for the unit string pattern: " [1/km^2]"
      cleaned_str <- sub(" \\[1/km\\^2\\]", "", as.character(val), fixed = TRUE)

      # Convert the resulting clean string to numeric
      return(as.numeric(cleaned_str))
    })

    # 2. Create tibble for this country
    tibble(
      Value = cleaned_values,
      Site = site_names
    )
  }, .id = "Country") %>%
    # Ensure the Country column is a factor for proper plotting order
    mutate(Country = factor(Country, levels = country_names_from_list)) %>%
    # Filter out rows where Value is NA (due to NULL input or failed string-to-numeric conversion)
    filter(!is.na(Value))

  # Print the transformed data structure for review
  cat("--- Transformed Data Frame Structure (NA/NULL filtered) ---\n")
  print(head(df_occurrences))
  cat("Total sites visualized:", nrow(df_occurrences), "\n\n")

  # --- 3. DATA VISUALIZATION (Grouped Violin Plot with Jittered Points) ---

  # Reorder countries by median value for better readability
  country_order <- df_occurrences %>%
    group_by(Country) %>%
    # FIX: Add na.rm = TRUE to explicitly handle missing values during aggregation
    summarise(MedianValue = median(Value, na.rm = TRUE)) %>%
    arrange(desc(MedianValue)) %>%
    pull(Country)

  df_occurrences$Country <- factor(df_occurrences$Country, levels = country_order)

  # Create the final ggplot visualization
  occurrence_plot <- ggplot(df_occurrences, aes(x = Country, y = Value)) +

    # 1. Violin Plot: Shows the overall distribution shape for each country
    geom_violin(
      aes(fill = Country),
      alpha = 0.6,
      show.legend = FALSE,
      width = 1.2
    ) +

    # 2. Jittered Points: Shows every individual site occurrence value
    geom_jitter(
      color = "black",
      size = 1.5,
      alpha = 0.7,
      width = 0.25 # Control the spread of the points
    ) +

    # 3. Median/Boxplot summary (Now uses pointrange and returns required aesthetics)
    stat_summary(
      fun.data = function(x) {
        # Calculate median
        y <- median(x, na.rm = TRUE)
        # Define ymin and ymax close to y to draw a strong median marker.
        # This satisfies the ymin/ymax aesthetic requirement.
        # A small fixed buffer (e.g., 0.1) is added for visual effect.
        data.frame(y = y, ymin = y - 0.1, ymax = y + 0.1)
      },
      geom = "pointrange", # Changed to pointrange (better than crossbar for a single marker)
      width = 0.3,         # Controls the horizontal width of the line segment
      size = 1,            # Controls the size of the point/line
      color = "darkred"    # Highlight the median value
    ) +

    # Labels and theming
    labs(
      title = "Site Occurrence Density Distribution by Country (Units Removed)",
      subtitle = paste("Countries (n=", length(country_order), ") ordered by median occurrence density | Total Sites:", nrow(df_occurrences)),
      y = "Median Density (Occurrences per km^2)",
      x = NULL
    ) +

    # IMPROVED: Logarithmic Y-axis to handle outliers and show distribution details
    # Now uses labels = scales::label_number() to display the untransformed values on the axis
    scale_y_log10(
      labels = scales::label_number()
    ) +
    theme_minimal(base_size = 14) +
    theme(
      # Rotate country labels for better fit on the x-axis
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      panel.grid.major.x = element_blank(), # Remove vertical grid lines
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "gray40")
    )

  # Display the plot
  return(occurrence_plot)

}
