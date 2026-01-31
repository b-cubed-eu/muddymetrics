# Technology Stack: muddymetrics

## Core Technologies
- **Programming Language:** R
- **Version Control:** Git & GitHub

## Libraries & Frameworks
- **Indicators:** `b3gbi` (Primary tool for calculating biodiversity metrics from GBIF data cubes)
- **Visualization:** `ggplot2` (Standardized maps and time-series plots using `viridis` scales)
- **Spatial Data:** `sf`, `terra` (Handling of WKT boundaries and spatial intersections)
- **Documentation:** `roxygen2` (In-code documentation for modular functions)
- **Testing:** `testthat` (Unit testing for indicator consistency)

## Data Architecture
- **Input:** GBIF Data Cubes (CSV/Parquet format), Ramsar Site Boundaries (WKT)
- **Intermediate/Output:** `.rds` (R-specific serialized data), `.png` (Visualizations)
- **Storage Strategy:** Strict separation of raw data (`inst/extdata`) and generated outputs (`output/`)
