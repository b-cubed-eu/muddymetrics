# Track Plan: Advanced Ecological Analysis

## Phase 1: Diversity and Rarity Metrics
- [ ] Task: Implement `calculate_alpha_diversity()` (Hill numbers, Shannon, Simpson) in `R/ecological_functions.R`.
- [ ] Task: Implement `calculate_evenness()` (Pielou, Williams).
- [ ] Task: Implement `calculate_rarity_metrics()` (Area/Abundance based).
- [ ] Task: Conductor - User Manual Verification 'Diversity and Rarity' (Protocol in workflow.md)

## Phase 2: Specialized Indicators Integration
- [ ] Task: Install and configure `pdindicatoR` and `impIndicatoR` packages.
- [ ] Task: Implement wrapper function `calculate_phylogenetic_diversity()` using `pdindicatoR`.
- [ ] Task: Implement wrapper function `calculate_invasive_impact()` using `impIndicatoR` and EICAT data.
- [ ] Task: Conductor - User Manual Verification 'Specialized Integration' (Protocol in workflow.md)

## Phase 3: European Regional Analysis
- [ ] Task: Create a workflow script `scripts/run_european_deep_dive.R` that filters for European sites and runs the advanced suite.
- [ ] Task: Generate comparative plots (e.g., Violin plots comparing diversity across European countries).
- [ ] Task: Conductor - User Manual Verification 'European Analysis' (Protocol in workflow.md)
