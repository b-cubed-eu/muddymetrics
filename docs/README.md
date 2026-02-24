# Deploying to GitHub Pages

This directory contains a static HTML gallery for browsing pre-generated plots.

## Quick Deploy

1. **Enable GitHub Pages:**
   - Go to your repository on GitHub
   - Settings → Pages
   - Source: Deploy from a branch
   - Branch: `main` (or your default branch)
   - Folder: `docs`
   - Click Save

2. **Wait 1-2 minutes** for deployment

3. Your site will be live at: `https://<your-username>.github.io/<repo-name>/`

## Rebuilding the Gallery

If you add new plots, regenerate the gallery:

```r
Rscript scripts/generate_static_gallery.R
```

This will update all HTML files in the `docs/` folder.

## Local Preview

To preview locally before pushing:

```r
# Option 1: Simple HTTP server
cd docs
python -m http.server 8000

# Option 2: Using R
library(servr)
httd("docs")
```

Then open http://localhost:8000 in your browser.

## Structure

```
docs/
├── index.html              # Home page - continents
├── africa/
│   └── index.html         # Countries in Africa
│   └── Country/
│       └── index.html     # Sites in Country
│           └── site_id/
│               └── index.html  # Plots for site
├── europe/
│   └── ...
└── ...
```

## Customization

Edit `scripts/generate_static_gallery.R` to customize:
- Colors (search for `#e94560`)
- Layout (grid columns, card styles)
- Navigation
