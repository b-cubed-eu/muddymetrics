# Track Specification: Advanced Ecological Analysis (European Deep Dive)

## Goal
Perform a deep-dive ecological analysis on European Ramsar sites (and other selected data-rich sites) using advanced biodiversity metrics. This extends the basic sufficiency analysis to include community composition and threat indicators.

## Core Indicators
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
