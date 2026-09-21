# muddymetrics

[![repo
status](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#wip)
[![Release](https://img.shields.io/github/release/b-cubed-eu/muddymetrics.svg?include_prereleases)](https://github.com/b-cubed-eu/muddymetrics/releases)
[![R-CMD-check](https://github.com/b-cubed-eu/muddymetrics/actions/workflows/R-CMD-check.yaml/badge.svg?branch=main)](https://github.com/b-cubed-eu/muddymetrics/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/b-cubed-eu/muddymetrics/branch/main/graph/badge.svg)](https://app.codecov.io/gh/b-cubed-eu/muddymetrics/)

Analysis code for the Ramsar case study in B3 deliverable D6.1, *Biodiversity change*. The
analysis tests whether open GBIF occurrence data are sufficient to monitor biodiversity trends
at Wetlands of International Importance, computes b3gbi temporal indicators for the sites that
pass, and illustrates invasive-species impact analysis at three well-sampled sites.

## Repository contents

- `scripts/`: the analysis workflow. See [`scripts/README.md`](scripts/README.md) for run order.
- `R/`, `man/`, `tests/`: helper functions packaged as `muddymetrics`. The analysis
  scripts are standalone and do not require the package to be installed.
- `docs/`: static gallery of per-site plots, published with GitHub Pages.

## Requirements

R with the following packages:

```r
install.packages(c("dplyr", "tidyr", "readr", "stringr", "stringi", "purrr", "data.table",
                   "sf", "units", "vegan", "callr", "withr",
                   "ggplot2", "scales", "patchwork", "gridExtra", "ggrepel", "svglite"))
remotes::install_github("b-cubed-eu/b3gbi")
remotes::install_github("b-cubed-eu/impIndicator")
```

## Data

Input data are not included in this repository:

- GBIF occurrence cubes at 100 m MGRS resolution, downloaded per continent from the GBIF web interface
- Ramsar site boundaries from the Ramsar Sites Information Service
- GRIIS national checklists of introduced and invasive species
- GIDIAS invasive-species impact records

The scripts expect these under `inst/extdata/`. See `scripts/README.md` for details.

## License

MIT. See [`LICENSE.md`](LICENSE.md). Citation details are in [`CITATION.cff`](CITATION.cff).

## Funding

This project receives funding from the European Union's Horizon Europe Research and Innovation
Programme (ID No 101059592).
