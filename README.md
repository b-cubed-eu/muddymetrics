# muddymetrics: Ramsar Biodiversity Indicator Pipeline

Evaluate Ramsar policy monitoring using open GBIF data and provide a global dashboard of biodiversity trends and indicators (richness, occupancy, evenness, rarity) for researchers and policymakers.

## Project Structure

- `R/`: Modular function definitions (Download, Indicator, Visualization, Utils).
- `scripts/`: Production scripts for data acquisition and heavy processing.
- `inst/extdata/`: Raw data, including Ramsar boundaries and GBIF cubes.
- `output/`: Generated results, including plots (.png) and processed data (.rds/.RData).
- `tests/`: Unit tests using the `testthat` framework.
- `main.R`: Master workflow script to orchestrate the pipeline.

## Data Ingestion Workflows

This project supports two primary ways to obtain GBIF data:

### 1. The "Continental" Workflow (Recommended for Global Analysis)
Use this for processing thousands of sites efficiently.
1.  **Download:** Run `scripts/download_continental_cubes.R`. This requests a massive 100m MGRS cube for an entire continent via the GBIF SQL API.
    *   **Resumable:** If the download (which can be >100GB) is interrupted, simply restart the script. It uses system `curl` to resume the transfer exactly where it left off.
    *   **Requirements:** ~200GB free disk space for ASIA or EUROPE.
2.  **Cut:** Run `scripts/split_data_cubes.R` to intersect the continental cube with your Ramsar polygons locally.
3.  **Analyze:** Use `calc_ramsar_indicator()` in `main.R`.

### 2. The "Site-by-Site" Workflow (Fallback / Single Site)
Use this for targeted updates or testing.
*   **Scripts:** `scripts/download_gbif_ramsar_by_site.R`.
*   **How:** Requests small, site-specific data cuts directly from GBIF.

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