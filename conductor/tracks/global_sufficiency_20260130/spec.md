# Track Specification: Global Data Sufficiency Assessment

## Goal
Implement a robust pipeline to calculate "Phase 1" data sufficiency indicators for all ~2,500 global Ramsar sites using GBIF data cubes. This establishes the baseline for how well Ramsar sites are represented in open biodiversity data.

## Core Indicators
1.  **Observed Richness:** Total number of unique species recorded.
2.  **Total Occurrences:** Raw count of occurrence records.
3.  **Cumulative Richness:** Species accumulation over time.
4.  **Density of Occurrences:** Occurrences per unit area (km²). Target threshold: **≥ 0.25 records/km²** (Troia & McManamay 2016).
5.  **Mean Year of Occurrence:** Average year of observation (proxy for data recency).
6.  **Inventory Completeness (Chao2):** Ratio of observed to estimated richness. Target threshold: **≥ 0.7**.
7.  **Survey Saturation (SAC Slope):** Slope of the Species Accumulation Curve. Target threshold: **≤ 0.10**.

## Requirements
- **Input:** Global GBIF occurrence cubes (CSV/Parquet) and Ramsar site boundaries (WKT).
- **Processing:** Iterate through all global sites, calculating the core indicators and the Troia 2016 "Moderate Threshold" suite.
- **Filtering:** Sites must be classified as "Data-Rich" (pass thresholds) or "Data-Poor" (quantify the monitoring gap).
- **Output:**
    - Standardized `.rds` files for each site/indicator.
    - Aggregated summary table (CSV) for global analysis.
    - Static plots for each indicator per site (time series and maps).
- **Performance:** Pipeline must handle ~2,500 sites efficiently (parallel processing recommended).

## Success Criteria
- Successful generation of all 5 indicators for >95% of valid Ramsar sites.
- Aggregated global summary table produced.
- Visualization generated for a representative sample of sites.
