# Product Guidelines: muddymetrics

## Development Philosophy
- **Reproducibility & Rigor:** Emphasis on reproducible scripts, detailed comments, and clear documentation for scientific verification.
- **B-Cubed Standards:** Adhere to the software development guides at [docs.b-cubed.eu](https://docs.b-cubed.eu/guides/software-development/).

## Visualization Standards
- **Standardized Color Schemes:** Use color-blind friendly and scientifically accurate scales (e.g., `viridis`) across all outputs.
- **Uniform Layouts:** Maintain consistent aspect ratios, labels, and legends to ensure comparability between different sites and regions.

## Data Management
- **Tidy Principles:** Strictly adhere to "tidy data" principles for all intermediate and final results.
- **Strict Separation:** Maintain clear separation between raw data (`inst/extdata`), intermediate processing, and final outputs (`output/`).
- **Versioning:** All code and relevant non-binary artifacts must be versioned and managed via GitHub.

## Collaborative Workflow
- **Modular Design:** Prioritize writing modular, reusable functions instead of monolithic scripts.
- **Unit Testing:** Implement unit tests for core indicator calculation logic to ensure consistency.
- **Transparent Documentation:** Use `roxygen2` for function documentation to ensure clear API definitions.
