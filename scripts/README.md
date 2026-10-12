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

1. Per-site cubes: `split_data_cubes_targeted.R`, `split_data_cubes_remaining.R`; then the 2024-2026
   update and one taxonomy for all years: `patch_site_cubes_update.R`, `harmonise_cube_taxonomy.R`
   (see "Updating the site cubes" below); then `check_site_cubes.R` before any analysis
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
onwards (`--last-year 2023` drops later years). Thresholds are applied in a separate step, so they can be varied
without recomputing.

The taxon subsets are birds, mammals, herpetofauna, insects, plants and fish. All taxa together
("full") is no longer screened by default, since the 2024-2026 update covers only these groups;
`--subsets full,aves,...` still includes it.

    Rscript scripts/run_sufficiency.R --cores 8            # writes output/sufficiency_v1/sufficiency_raw.csv
    Rscript scripts/run_sufficiency.R --cores 8 --grid 1km # the same on 1 km cells, in output/sufficiency_v1_1km/
    Rscript scripts/apply_sufficiency_thresholds.R         # writes sufficiency_results.csv, sites_passing.csv
    Rscript scripts/apply_sufficiency_thresholds.R --precision-rule median   # alternative precision rule
    Rscript scripts/apply_sufficiency_thresholds.R --in output/sufficiency_v1_1km/sufficiency_raw.csv
    Rscript scripts/apply_sufficiency_thresholds.R --coverage-rule informative --min-informative-years 10
    Rscript scripts/check_sufficiency_against_d61.R        # optional comparison with the D6.1 results
    Rscript -e 'testthat::test_file("scripts/tests/test_sufficiency_criteria.R")'

Every run also reports, per site x taxon, the share of rows and of occurrences whose uncertainty circle is
larger than 0.1 x the site area (`share_rows_vague`, `share_occ_vague`; among rows with an uncertainty value).
`--drop-vague` leaves those records out before any criterion is computed (sensitivity check for the precision
criterion; use its own `--out` folder).

The runner is resumable and runs each site in its own R process. All criterion functions are in
`scripts/lib/sufficiency_criteria.R`. Requires data.table, dplyr, sf, iNEXT and testthat.

| Criterion | Definition | Pass |
|---|---|---|
| Density | occurrences / geodesic site area (km²) | ≥ 0.25 |
| Chao2 completeness | observed / bias-corrected Chao2 richness; species x cell incidence (100 m, or 1 km with `--grid 1km`), years pooled | ≥ 0.70 |
| Coverage deficit | 1 − mean annual incidence sample coverage (iNEXT), over years with at least 2 occupied cells and more than 5 species x cell entries (`--coverage-rule informative`: over informative years only, see below) | ≤ 0.10 |
| Temporal decoupling | Cha (2007) distance between min-max scaled annual total occurrences and observed richness | ≥ elbow threshold |
| Spatial precision | occurrence-weighted mean uncertainty area (π r²) / site area (`--precision-rule median`: weighted median; `d61`: the D6.1 site-level definition) | ≤ 0.10 |

Only Chao2 completeness and annual coverage use grid cells as sampling units, so only they change
with `--grid`.

An informative year has at least 20 occupied cells and 50 species x cell incidences
(`INFORMATIVE_MIN_CELLS`, `INFORMATIVE_MIN_INCIDENCES` in `scripts/lib/sufficiency_criteria.R`). The limits come
from `scripts/check_informative_years.R`, which subsamples well-sampled 1 km site-years: in years below
them the estimated annual coverage is too often 0.90 or more while the true coverage is below 0.85.

A criterion that cannot be computed (for example Chao2 with fewer than two occupied cells) counts as a
fail. Species are identified by scientific name, because the cubes carry two GBIF key systems and the
same species can appear under two keys. Herpetofauna includes Amphibia, Reptilia, Squamata, Testudines,
Crocodylia and Rhynchocephalia. Fish are chordates in the fish classes (`FISH_CLASSES`) plus chordates
without a class, which in the old GBIF backbone are the ray-finned fishes.

The threshold step also prints how many combinations pass under each coverage and precision rule.

## Indicators, version 1

    Rscript scripts/run_indicators.R --cores 8                                    # point estimates, all evaluated combinations
    Rscript scripts/run_indicators.R --cores 8 --mode ci --nboot 100 --grids 1km  # bootstrap intervals, passing combinations
    Rscript scripts/summarise_indicators.R --grid 1km --level 0.9 --ci hill_ci_1km_b100.csv   # trends; repeat per grid and level
    Rscript scripts/check_multiverse.R --precision-rule median     # verdicts across analysis choices
    Rscript scripts/check_detectability.R --precision-rule median  # smallest detectable trend

For every site x taxon combination the screen evaluated, `worker_indicators.R` computes total
occurrences, observed and cumulative richness (by scientific name), occupied cells, annual sample
coverage and the Hill numbers q = 0, 1, 2 at 90% and 95% sample coverage, on the native 100 m MGRS
grid and on 1 km MGRS cells, plus the MAGV reliability series (`scripts/lib/magv.R`). Combinations
that pass all five criteria under any precision rule also get bootstrap intervals for the Hill numbers.

The Hill numbers (`scripts/lib/indicators_v1.R`) use b3gbi's coverage-based estimator (iNEXT
`estimateD` through b3gbi's wrapper) on yearly incidence frequencies, which gives the same values as
`b3gbi::hill0_ts()` / `hill1_ts()` / `hill2_ts()` on the native grid without building the species x
cell matrix. Each estimate records the sample size needed to reach the coverage level, the observed
number of cells, iNEXT's method (rarefaction, observed or extrapolation) and whether its interval
could be computed. 1 km cells are the MGRS 1 km squares containing the 100 m cells.

By default `summarise_indicators.R` keeps Hill estimates only from informative years and only where reaching the
coverage level needs at most twice the observed number of cells (`--hill-years all --max-extrapolation Inf`
keeps all), and `--to-year 2023` ends the series in 2023. `--in` selects the indicator folder.

`run_indicators.R` takes `--screen`, `--coverage-rule` and `--out` so the indicators can follow another
screen run. Two options add to an earlier interval run: `--combos` (a CSV of site_id and subset)
computes intervals for particular combinations, and `--years 2024-2025` for particular years only.
Each year is estimated on its own, with its own random seed, so intervals for added years are the
same as in a run over all years.

`summarise_indicators.R` fits each annual series by ordinary least squares and by generalised
least squares with AR(1) errors (REML, `nlme::gls`), and, with `--ci`, both again with each year
weighted by its bootstrap variance. A trend counts as up or down when p ≤ `--alpha` (default 0.05).
`check_multiverse.R` compares the Hill-number verdicts across Hill orders, grids, coverage levels and
trend models; `check_detectability.R` gives the linear change per decade that a trend test over a
chosen window would detect with the given power. Both also need data.table, and the trend models nlme.

## Updating the site cubes (2024-2026 data, one taxonomy)

The June 2026 update (GBIF download 0074420-260519110011954, new GBIF backbone) did not reach
every site cube: sites in countries with folders under more than one continent, sites re-split on
16-17 June and sites with a dot in their name were skipped, and reptiles were dropped everywhere.
Two scripts repair this; both copy every file they change to a backup folder first.

1. `patch_site_cubes_update.R` rebuilds the rows from 2024 on from the update file, assigning cells
   to sites by intersecting them with the site outlines, as in the original split. Rows before 2024
   are kept. Several update files can be given comma-separated, each with the taxa taken from it
   (`--update-groups`): the main update for all its groups, and the fish from a second download (GBIF
   download 0014182-260928105237408, which holds all bony vertebrates). `--to-year 2025` leaves out
   2026, a partial year. Check first, then run:
   ```
   Rscript scripts/patch_site_cubes_update.R --sites site_67,site_838 --out-dir _local_data/cube_patch_check \
     --update <main>,<fish> --update-groups all,fish --to-year 2025
   Rscript scripts/patch_site_cubes_update.R --all --update <main>,<fish> --update-groups all,fish --to-year 2025
   ```
   or, to check the result before any cube is changed, with the updates read only once:
   ```
   Rscript scripts/patch_site_cubes_update.R --all --out-dir _local_data/cube_patch_all --update <main>,<fish> \
     --update-groups all,fish --to-year 2025
   Rscript scripts/patch_site_cubes_update.R --install-from _local_data/cube_patch_all --dry-run
   Rscript scripts/patch_site_cubes_update.R --install-from _local_data/cube_patch_all
   ```
   The install step copies each original to the backup folder, verifies the copies, and writes
   `install_log_<date>_<time>.csv` in the work folder.
   For cubes that already had the update, the summary reports whether the rows rebuilt from the
   main update are identical apart from reptiles, and how many rows came from the other files.
2. `harmonise_cube_taxonomy.R` gives the older rows (old backbone, numeric keys) the names and keys
   of the new backbone, so that a species keeps one name across years. It builds a name table (from
   the update rows, then the GBIF species-match API; needs internet, resumable) and applies it.
   Only names in the screened taxon groups are matched (`--groups all` matches every name):
   ```
   Rscript scripts/harmonise_cube_taxonomy.R --build --limit 200
   Rscript scripts/harmonise_cube_taxonomy.R --build
   Rscript scripts/harmonise_cube_taxonomy.R --apply
   ```

Functions: `lib/cube_update.R`, `lib/taxonomy_match.R`; tests: `tests/test_cube_update.R`.

## Checking the site cubes before analysis

`check_site_cubes.R` rebuilds, from the sources alone, what every site cube should contain, and
compares it row by row with the cubes on disk. It never changes a cube. Its code is separate from the
code that built the cubes, so that an error there is not repeated in the check:

- outlines come from the Ramsar shapefile; the .wkt files are compared with it separately;
- the continental cubes give the years before 2024 (all of them, for every site), the update files
  2024 on;
- every line of every source is accounted for, and every valid cell is tested against every outline
  (no pre-selection by 100 km square).

A cell belongs to a site when the south-west corner of the 100 m cell (the point
`mgrs::mgrs_to_latlng()` returns) lies in or on the outline, the rule the cubes were built with.

```
Rscript scripts/check_site_cubes.R --step inventory                                  # files and outlines (minutes)
Rscript scripts/check_site_cubes.R --step reference --update <main>,<fish> --update-groups all,fish --to-year 2025 --cores 8
                                                     # read all sources (hours; resumable)
Rscript scripts/check_site_cubes.R --step compare --tag after_patch --to-year 2025
Rscript scripts/check_site_cubes.R --step report --tag after_patch
Rscript scripts/check_site_cubes.R --step compare --tag after_taxonomy --to-year 2025 --name-table _local_data/taxonomy/name_table.csv
Rscript scripts/check_site_cubes.R --step report --tag after_taxonomy
```

The reference step is the slow one and is done once; the compare step can be repeated (for example
after the patch and again after the taxonomy harmonisation, when the reference is renamed with the
name table before comparing). Output goes to `_local_data/cube_check/`; `REPORT_<tag>.md` explains
each check and ends with a verdict. Functions: `lib/cube_check.R`, `lib/cube_check_report.R`; tests
(with a made-up world containing known errors): `tests/test_check_cubes.R`.
Rows that match their source row on every column except the higher-taxon counts (kingdomcount to
genuscount) are counted separately, not as differences: the screen and the indicators do not use
these columns.

The check of October 2026 found cubes whose rows before 2024 are incomplete (cubes split on 16-17
June 2026 that stopped reading early or used another continent's cube, sites whose outlines are
invalid in planar geometry, rows filed under a neighbouring continent). `rebuild_site_cubes_from_reference.R`
replaces the rows before 2024 of these cubes with the rows the reference step collected, keeps the
rows from 2024 on, and writes new cubes for sites with data but no cube. Cubes filed under another
continent than the one whose continental cube holds most of their older rows (such as the Corsican
Etang de Biguglia under `ramsar_site_data_100m_southamerica/France`) are moved to that continent's
folder; the v1 scripts find cubes by site id (`discover_site_cubes()`), so only the continent label
in their output changes. The rebuilt cubes go to a
separate folder and are put into place with the patch script's install step; the check is then run
again. For the rebuilt sites the check is no longer independent of the build, since both use the
reference rows.
```
Rscript scripts/rebuild_site_cubes_from_reference.R --out-dir _local_data/cube_rebuild --cores 6
Rscript scripts/patch_site_cubes_update.R --install-from _local_data/cube_rebuild \
  --backup-dir _local_data/cube_backup_before_rebuild --work-dir _local_data/cube_rebuild_work
Rscript scripts/check_site_cubes.R --step compare --tag after_rebuild --update <main>,<fish> --update-groups all,fish --to-year 2025
Rscript scripts/check_site_cubes.R --step report --tag after_rebuild
```
Functions: `lib/cube_rebuild.R`.

