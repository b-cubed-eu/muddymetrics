# Track Plan: Global Data Sufficiency Assessment

## Phase 1: Pipeline Setup and Data Prep
- [x] Task: Preprocess global Ramsar shapefile and generate WKT boundaries.
- [x] Task: Organize site-level data directories.
- [ ] Task: Conductor - User Manual Verification 'Pipeline Setup' (Protocol in workflow.md)

## Phase 2: Indicator Implementation & Validation

- [x] Task: Calculate Observed Richness for global sites.

- [x] Task: Calculate Total Occurrences for global sites.

- [x] Task: Calculate Cumulative Richness for global sites.

- [ ] Task: Validate existing `overall_density.RData` files and ensure consistency with the **0.25 records/km²** threshold.

- [ ] Task: Implement "Mean Year of Occurrence" across all sites to assess data recency.

- [ ] Task: Implement **Chao2** and **SAC Slope** calculation functions in `R/`.

- [ ] Task: Conductor - User Manual Verification 'Indicator Implementation' (Protocol in workflow.md)



## Phase 3: Global Execution and Aggregation

- [ ] Task: Create a script `R/aggregate_global_results.R` to consolidate all site-level `.RData` files into a single master dataset.

- [ ] Task: Generate a global summary CSV `output/global_sufficiency_summary.csv` with key metrics for all ~2,500 sites.

- [ ] Task: Perform "Data Gap Analysis" and classify sites into **Data-Rich** vs **Data-Poor** cohorts based on the Troia 2016 criteria (Density ≥ 0.25, Chao2 ≥ 0.7, Slope ≤ 0.1).

- [ ] Task: Conductor - User Manual Verification 'Global Execution' (Protocol in workflow.md)
