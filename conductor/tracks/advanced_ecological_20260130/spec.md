# Track Specification: Advanced Ecological Analysis (European Deep Dive)

## Goal
Perform a deep-dive ecological analysis on European Ramsar sites and any global sites that meet the **Troia 2016 "Moderate Threshold"** for sampling completeness. This ensures that advanced metrics like evenness and turnover are only calculated for sites where the underlying data is scientifically robust.

## Site Cohort Definition
- **Geographic:** All European sites (for comparative regional analysis).
- **Sufficiency:** Any global site passing the "Data-Rich" filter from Track 2:
    - Density ≥ 0.25 records/km².
    - Chao2 Inventory Completeness ≥ 0.7.
    - SAC Slope ≤ 0.10.
1.  **Diversity & Evenness:** Pielou’s Evenness, Williams’ Evenness, Hill-Shannon Diversity, Hill-Simpson Diversity.
2.  **Rarity:** Area-Based Rarity, Abundance-Based Rarity.
3.  **Turnover:** Occupancy Turnover over time.
4.  **Specialized:**
    - **Phylogenetic Diversity:** Using `pdindicatoR`.
    - **Invasive Species Impact:** Using `impIndicatoR` (EICAT framework).

## Requirements
- **Scope:** Restricted to European sites (Continental Analysis) and selected high-quality global sites (Deep Dive).
- **Dependencies:** `b3gbi`, `pdindicatoR`, `impIndicatoR`, `vegan` (likely dependency for diversity metrics).
- **Output:** Detailed ecological profiles for selected sites.

## Success Criteria
- Advanced indicators calculated for the defined subset of sites.
- Successful integration of `pdindicatoR` and `impIndicatoR` packages.
- Comparative analysis of European sites completed.
