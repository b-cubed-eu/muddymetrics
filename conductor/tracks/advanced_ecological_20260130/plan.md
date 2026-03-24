# Track Plan: Advanced Ecological Analysis

## Phase 1: Diversity and Rarity Metrics
- [x] Task: Implement `calculate_alpha_diversity()` (Hill numbers, Shannon, Simpson) in `R/ecological_functions.R`.
  - Created: `R/ecological_functions.R` with Hill 0/1/2, Shannon, Simpson metrics.
- [x] Task: Implement `calculate_evenness()` (Pielou, Williams).
  - Added: Pielou's J and Williams' W calculations.
- [x] Task: Implement `calculate_rarity_metrics()` (Area/Abundance based).
  - Added: Area-based rarity and abundance-based rarity metrics.
- [ ] Task: Conductor - User Manual Verification 'Diversity and Rarity' (Protocol in workflow.md)

## Phase 2: Specialized Indicators Integration
- [x] Task: Install and configure `pdindicatoR` and `impIndicatoR` packages.
  - Installed from https://b-cubed-eu.r-universe.dev
- [x] Task: Implement wrapper function `calculate_phylogenetic_diversity()` using `pdindicatoR`.
  - Created: R/specialized_indicators.R with pdindicatoR wrapper
- [x] Task: Implement wrapper function `calculate_invasive_impact()` using `impIndicator` and EICAT data.
  - Added: impIndicator wrapper with EICAT integration
- [ ] Task: Conductor - User Manual Verification 'Specialized Integration' (Protocol in workflow.md)

## Phase 3: European Regional Analysis
- [x] Task: Create a workflow script `scripts/run_european_deep_dive.R` that filters for European sites and runs the advanced suite.
  - Created: Complete workflow with alpha diversity, inventory completeness, mean year calculations
- [x] Task: Generate comparative plots (e.g., Violin plots comparing diversity across European countries).
  - Added: plot_diversity_by_country() and plot_completeness_by_country() functions
- [ ] Task: Conductor - User Manual Verification 'European Analysis' (Protocol in workflow.md)
