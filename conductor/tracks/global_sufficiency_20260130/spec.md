# Track Specification: Global Data Sufficiency Assessment

## Goal
Implement a robust pipeline to calculate "Phase 1" data sufficiency indicators for all ~2,500 global Ramsar sites using GBIF data cubes. This establishes the baseline for how well Ramsar sites are represented in open biodiversity data.

## Core Indicators
1.  **Observed Richness:** Total number of unique species recorded.
2.  **Total Occurrences:** Raw count of occurrence records.
3.  **Cumulative Richness:** Species accumulation over time.
4.  **Density of Occurrences:** Occurrences per unit area (km²).
5.  **Mean Year of Occurrence:** Average year of observation (proxy for data recency).

## Requirements
- **Input:** Global GBIF occurrence cubes (CSV/Parquet) and Ramsar site boundaries (WKT).
- **Processing:** Iterate through all global sites, calculating the 5 core indicators using `b3gbi`.
- **Output:**
    - Standardized `.rds` files for each site/indicator.
    - Aggregated summary table (CSV) for global analysis.
    - Static plots for each indicator per site (time series and maps).
- **Performance:** Pipeline must handle ~2,500 sites efficiently (parallel processing recommended).

## Success Criteria
- Successful generation of all 5 indicators for >95% of valid Ramsar sites.
- Aggregated global summary table produced.
- Visualization generated for a representative sample of sites.
