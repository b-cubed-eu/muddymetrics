# muddymetrics: Ramsar Biodiversity Indicator Pipeline

Evaluate Ramsar policy monitoring using open GBIF data and provide a global dashboard of biodiversity trends and indicators (richness, occupancy, evenness, rarity) for researchers and policymakers.

## Project Structure

- `R/`: Modular function definitions (Download, Indicator, Visualization).
- `scripts/`: Long-running scripts for data acquisition and heavy processing.
- `inst/extdata/`: Raw data, including Ramsar boundaries and GBIF cubes.
- `output/`: Generated results, including plots (.png) and processed data (.rds/.RData).
- `tests/`: Unit tests using the `testthat` framework.
- `main.R`: Master workflow script to orchestrate the pipeline.

## Getting Started

### Prerequisites

- R (version 4.3.1 recommended)
- `b3gbi` package (installed from GitHub: `b-cubed-eu/b3gbi`)
- GBIF credentials (for data download)

### Usage

1. **Setup:** Ensure all dependencies are installed.
2. **Configuration:** Set your GBIF credentials in your `.Renviron` or via `options()`.
3. **Execution:** Run `main.R` to execute the example workflow.
4. **Testing:** Run `testthat::test_dir('tests/testthat/')` to verify the installation.

## Core Modules

- **Download:** `get_gbif_predicates()` for robust GBIF API queries.
- **Indicators:** `calculate_ramsar_metric()` wrapper for `b3gbi` workflows.
- **Visualization:** `save_ramsar_plot()` for standardized viridis-styled plotting.
- **Batch Processing:** `calc_ramsar_indicator()` for iterating over thousands of sites.

## Development

This project follows a spec-driven development framework (Conductor) with a focus on scientific rigor, TDD, and modularity.
