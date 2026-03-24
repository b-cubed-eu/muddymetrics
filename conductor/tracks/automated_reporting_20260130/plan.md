# Track Plan: Automated Report Generation (Dashboard)

## Phase 1: Automated Site Reports
- [x] Task: Design an RMarkdown template `inst/templates/site_report_template.Rmd`.
  - Created: Full template with site info, data sufficiency, plots
- [x] Task: Implement a "Report Orchestrator" script to loop through site folders and render individual PDFs.
  - Created: R/report_orchestrator.R with generate_site_reports() and generate_summary_report()
- [ ] Task: Conductor - User Manual Verification 'Site Reports' (Protocol in workflow.md)

## Phase 2: Global Dashboard Prototype
- [x] Task: Set up a basic Shiny application structure in `dashboard/`.
  - Created: dashboard/app.R with global.R
- [x] Task: Implement global/continental map view using the master summary dataset.
  - Added: leaflet map with site markers and filters
- [x] Task: Implement site-specific "Deep Dive" views within the dashboard.
  - Added: Site details tab with time series plots
- [ ] Task: Conductor - User Manual Verification 'Dashboard Prototype' (Protocol in workflow.md)

## Phase 3: Final Polish and Deployment
- [ ] Task: Apply B3 branding and CSS styling to all reports and the dashboard.
- [ ] Task: Optimize dashboard performance for large site lists (using `leaflet` or `mapdeck`).
- [ ] Task: Conductor - User Manual Verification 'Final Polish' (Protocol in workflow.md)
