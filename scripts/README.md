# Analysis scripts

Scripts for the Ramsar Wetlands of International Importance case study (B3 deliverable D6.1).

Run all scripts from the repository root, e.g. `Rscript scripts/run_subgroup_data_sufficiency_parallel.R`.
Runner scripts call their workers by the relative path `scripts/<worker>.R`.

## Input data

Input data are not included in this repository. Scripts expect GBIF occurrence cubes,
Ramsar site boundaries, GRIIS checklists and GIDIAS impact records under `inst/extdata/`
(GIDIAS at `inst/extdata/GIDIAS/GIDIAS_machine_read.csv`).

## Run order

1. Per-site cubes: `split_data_cubes_targeted.R`, `split_data_cubes_remaining.R`
2. Data-sufficiency screen: `run_subgroup_data_sufficiency_parallel.R` (calls `worker_subgroup.R`
   and `compute_thresholds.R`), then `global_precision_update.R` and `global_precision_update_corrected.R`
3. Summary tables: `generate_stats_tables.R`, `final_stats.R`
4. Indicators for passing sites: `run_full_analysis.R` (calls `worker_full_analysis.R`), then
   `make_part2_summaries.R`
5. Figures: `generate_report_figures.R`
6. Invasive species: `generate_invasive_vignette.R`, `generate_all_sites_invasive.R`,
   `generate_effort_standardised_invasive.R`, `generate_coverage_standardised_invasive.R`,
   `generate_composition_metric_invasive.R`

Outputs are written to `output/` and `report_figures/`.
