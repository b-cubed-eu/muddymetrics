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

## Data-sufficiency screen, version 1

A single worker computes all criteria for every site x taxon combination, from records dated 1990
onwards. Thresholds are applied in a separate step, so they can be varied without recomputing.

    Rscript scripts/run_sufficiency.R --cores 8            # writes output/sufficiency_v1/sufficiency_raw.csv
    Rscript scripts/apply_sufficiency_thresholds.R         # writes sufficiency_results.csv, sites_passing.csv
    Rscript scripts/apply_sufficiency_thresholds.R --precision-rule median   # alternative precision rule
    Rscript scripts/check_sufficiency_against_d61.R        # optional comparison with the D6.1 results
    Rscript -e 'testthat::test_file("scripts/tests/test_sufficiency_criteria.R")'

The runner is resumable and runs each site in its own R process. All criterion functions are in
`scripts/lib/sufficiency_criteria.R`. Requires data.table, dplyr, sf, iNEXT and testthat.

| Criterion | Definition | Pass |
|---|---|---|
| Density | occurrences / geodesic site area (km²) | ≥ 0.25 |
| Chao2 completeness | observed / bias-corrected Chao2 richness; species x 100 m cell incidence, years pooled | ≥ 0.70 |
| Coverage deficit | 1 − mean annual incidence sample coverage (iNEXT), over years with at least 2 occupied cells and more than 5 species x cell entries | ≤ 0.10 |
| Temporal decoupling | Cha (2007) distance between min-max scaled annual total occurrences and observed richness | ≥ elbow threshold |
| Spatial precision | occurrence-weighted mean uncertainty area (π r²) / site area (`--precision-rule median`: weighted median; `d61`: the D6.1 site-level definition) | ≤ 0.10 |

A criterion that cannot be computed (for example Chao2 with fewer than two occupied cells) counts as a
fail. Species are identified by scientific name, because the cubes carry two GBIF key systems and the
same species can appear under two keys. Herpetofauna includes Amphibia, Squamata, Testudines,
Crocodylia and Rhynchocephalia.

The threshold step also prints how many combinations pass under each coverage and precision rule.

## Indicators for passing combinations, version 1

    Rscript scripts/run_indicators.R --cores 2     # writes output/indicators_v1/
    Rscript scripts/summarise_indicators.R         # writes indicator_summary.csv, trend_direction.csv

`run_indicators.R` runs `worker_indicators.R` for every combination in
`output/sufficiency_v1/sites_passing.csv`: all b3gbi time-series indicators (normal-approximation
bootstrap intervals, 1000 replicates) and the MAGV reliability series (`scripts/lib/magv.R`). It
requires b3gbi 1.0.0 or later and harmonises species keys by scientific name before building the
cube, so a species recorded under two GBIF keys is counted once. Indicators that fail with intervals
are rerun without them, and the outcome of every indicator is listed in
`output/indicators_v1/status.csv`. Taxonomic distinctness queries GBIF through taxize, so it needs
internet access. Requires b3gbi, taxize, rnaturalearthdata, data.table and dplyr.
