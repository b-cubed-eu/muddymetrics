# muddymetrics: Ramsar Biodiversity Indicator Pipeline

[![repo
status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#wip)
[![Release](https://img.shields.io/github/release/b-cubed-eu/muddymetrics.svg?include_prereleases)](https://github.com/b-cubed-eu/muddymetrics/releases)
[![R-CMD-check](https://github.com/b-cubed-eu/muddymetrics/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/b-cubed-eu/muddymetrics/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/b-cubed-eu/muddymetrics/branch/main/graph/badge.svg)](https://app.codecov.io/gh/b-cubed-eu/muddymetrics/)


Evaluate Ramsar policy monitoring using open GBIF data and provide a global dashboard of biodiversity trends and indicators (richness, occupancy, evenness, rarity) for researchers and policymakers.

## Project Structure

- `R/`: Modular function definitions (Download, Indicator, Visualization, Utils).
- `scripts/`: Analysis workflow for the Ramsar case study. See `scripts/README.md` for run order.
- `inst/extdata/`: Raw data, including Ramsar boundaries and GBIF cubes.
- `output/`: Generated results, including plots (.png) and processed data (.rds/.RData).
- `tests/`: Unit tests using the `testthat` framework.

## Data

Continental GBIF occurrence cubes (100 m MGRS grid) were downloaded from the GBIF web interface and split into per-site cubes with `scripts/split_data_cubes_targeted.R` and `scripts/split_data_cubes_remaining.R`. Input data are not included in this repository; see `scripts/README.md` for the expected locations under `inst/extdata/`.

## Getting Started

### Prerequisites
- R (version 4.3.1 recommended)
- `b3gbi` package (installed from GitHub: `b-cubed-eu/b3gbi`)
- **System Curl:** Ensure `curl` is available in your system path (Standard on Windows 10+).

### Configuration
Set your GBIF credentials in your `.Renviron` file:
```R
GBIF_USER="your_username"
GBIF_PWD="your_password"
GBIF_EMAIL="your_email@example.com"
```

## Core Modules
- **Download:** `get_gbif_predicates()` for API queries.
- **Manager:** `download_robust()` handles resumable 100GB+ file transfers.
- **Indicators:** `calculate_ramsar_metric()` wrapper for `b3gbi`.
- **Batch:** `calc_ramsar_indicator()` handles iteration over thousands of files.

## Development
This project follows a spec-driven development framework (Conductor) with a focus on scientific rigor, TDD, and modularity. Run `testthat::test_dir('tests/testthat/')` to verify the installation.
