# muddymetrics 0.2.0

* Added the analysis workflow for the Ramsar case study (B3 deliverable D6.1) in `scripts/`,
  with run order and input data documented in `scripts/README.md`.
* `calculate_chao2()` now estimates Chao2 from species incidence across occupied grid cells,
  and returns `NA` when fewer than two cells are occupied.
* Removed superseded download and processing scripts and `main.R`.
* Rewrote the README to describe the analysis workflow, requirements and input data.
* Added `.gitattributes` to normalise line endings.

# muddymetrics 0.1.0

* Initial release.
* Modularized codebase into functional components in `R/`.
* Moved monolithic scripts to `scripts/`.
* Implemented automated batch indicator calculation with `calc_ramsar_indicator()`.
* Added robust, resumable continental GBIF download logic.
* Established comprehensive unit test suite with >80% coverage.
* Configured GitHub Actions for automated `R CMD check` and coverage reporting.
