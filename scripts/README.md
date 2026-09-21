# Analysis scripts

Scripts for the Ramsar Wetlands of International Importance case study (B3 deliverable D6.1).

Run all scripts from the repository root, e.g. `Rscript scripts/run_subgroup_data_sufficiency_parallel.R`.
Runner scripts call their workers by the relative path `scripts/<worker>.R`.

## Input data

Input data are not included in this repository. Continental GBIF occurrence cubes (100 m MGRS grid)
were downloaded from the GBIF web interface. Scripts expect these cubes,
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

`generate_static_gallery.R` builds the GitHub Pages gallery in `docs/` from per-site plots in `output/`.

## Revised temporal and density criteria

`run_temporal_only.R` (calling `worker_temporal_only.R`) recomputes the temporal-decoupling
criterion from annual rather than cumulative series and records per-combination occurrence
totals. `merge_and_regate_v2.R` then recomputes occurrence density from site area in km² and
re-applies the five criteria. Run after step 2:

    Rscript scripts/run_temporal_only.R
    Rscript scripts/merge_and_regate_v2.R
