# Plain-language report of scripts/check_site_cubes.R.

chk_write_report <- function(out, tag, from_year) {
  cmp <- file.path(out, paste0("compare_", tag))
  rd <- function(f) if (file.exists(f)) data.table::fread(f, colClasses = list(character = "site_id")) else NULL
  inv <- rd(file.path(out, "inventory.csv")); oc <- rd(file.path(out, "outline_check.csv"))
  src <- if (file.exists(file.path(out, "reference", "sources.csv"))) data.table::fread(file.path(out, "reference", "sources.csv")) else NULL
  st <- rd(file.path(cmp, "sites.csv")); yg <- rd(file.path(cmp, "year_group.csv.gz"))
  L <- character()
  add <- function(...) L <<- c(L, ...)
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  tab <- function(d) {
    d <- as.data.frame(d)
    for (j in seq_along(d)) d[[j]] <- if (is.numeric(d[[j]]) && names(d)[j] != "year") fmt(d[[j]]) else as.character(d[[j]])
    c(paste0("| ", paste(names(d), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
      apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  }
  verdict <- character()

  add(sprintf("# Site cube check (%s)", tag), "", sprintf("Written %s.", format(Sys.time(), "%d %B %Y, %H:%M")), "")
  add(paste0("This report compares every site cube with what its sources say it should contain. Sources: the Ramsar ",
             "shapefile for the outlines, the continental cubes for the years before ", from_year,
             ", and the update files for ", from_year, " on. A cell counts as part of a site when the south-west corner ",
             "of the 100 m cell lies in or on the site outline."), "")

  # 1. Files and outlines
  add("## 1. Files and outlines", "")
  if (!is.null(inv)) {
    n_shp <- sum(inv$in_shapefile)
    add(sprintf("The shapefile has %s sites. Cube files were found for %s of them and outline (.wkt) files for %s.",
                fmt(n_shp), fmt(sum(inv$in_shapefile & inv$n_cubes > 0)), fmt(sum(inv$in_shapefile & inv$n_wkt > 0))), "")
    iss <- inv[issues != ""]
    if (nrow(iss)) {
      cnt <- data.table::data.table(issue = unlist(strsplit(iss$issues, "; ")))[, .(sites = .N), by = issue][order(-sites)]
      add("Issues found (a site can have more than one):", "", tab(cnt), "",
          "All sites with an issue are listed in `inventory.csv`. A site without a cube file is only a problem if",
          "its sources hold data for it (see section 3).", "")
      if (any(!grepl("^no cube file$", iss$issues))) verdict <- c(verdict, "file issues other than missing cubes")
    } else add("No issues: every site has exactly one cube file and one outline file.", "")
    nts <- inv[notes != ""]
    if (nrow(nts)) {
      cnt <- data.table::data.table(note = unlist(strsplit(nts$notes, "; ")))[, .(sites = .N), by = note][order(-sites)]
      add("Notes (not counted as problems; listed in `inventory.csv`):", "", tab(cnt), "")
    }
  }
  if (!is.null(oc)) {
    cmpd <- oc[status == "compared"]
    add(sprintf("Outline files compared with the shapefile: %s (could not compare: %s).", fmt(nrow(cmpd)), fmt(sum(oc$status != "compared"))), "")
    q <- stats::quantile(cmpd$mean_shift_m, c(0.5, 0.9, 0.99, 1), na.rm = TRUE)
    add(sprintf("Average distance between the outline file's boundary and the shapefile's: median %.1f m, 90%% of sites below %.1f m, 99%% below %.1f m, largest %s m.",
                q[1], q[2], q[3], fmt(round(q[4]))),
        "A few metres means the same outline with rounded coordinates; a cell is 100 m wide.", "")
    big <- cmpd[mean_shift_m > 10][order(-mean_shift_m)]
    if (nrow(big)) {
      add(sprintf("%s sites have outline files whose boundary lies on average more than 10 m from the shapefile's (largest 25",
                  fmt(nrow(big))), "below; all in `outline_check.csv`). Their cubes were split with these outlines, so cells near the",
          "boundary can differ from this check, which uses the shapefile.", "",
          tab(utils::head(big[, .(site_id, mean_shift_m = round(mean_shift_m), share_of_area_differing = signif(rel_symdiff, 2),
                                  area_ratio = signif(area_ratio, 4))], 25)), "")
      verdict <- c(verdict, "outline files that differ from the shapefile")
    }
  }

  # 2. Sources
  add("## 2. Source files", "")
  if (!is.null(src)) {
    s <- src[, .(file = basename(source), role = ifelse(role == "before", paste("years before", from_year), paste(from_year, "on")),
                 taxa = groups, lines, read_completely, wrong_fields = n_wrong_fields, empty = n_empty,
                 invalid_cell = n_invalid_code, no_coordinates = n_no_coords, rows_in_sites = rows_kept,
                 rows_other_years = rows_other_period, rows_other_taxa, hours = round(minutes / 60, 1))]
    add(tab(s), "",
        "`read_completely` says the byte ranges were read without gaps or overlaps, to the end of the file.",
        "`wrong_fields` are lines that are not 25 tab-separated fields; they cannot be placed and would be missing",
        "from any cube. `invalid_cell` and `no_coordinates` count cell codes (once per block of lines read) that are not valid 100 m MGRS",
        "codes or could not be converted. `rows_other_years` are rows in a site from the other period (continental",
        sprintf("rows from %d on, update rows outside the update years); they are left out of the cubes by design.", from_year),
        "`rows_other_taxa` are rows of an update file outside its `taxa`; the classes left out are listed below, so",
        "that a class missing from a taxon definition shows up here.", "")
    oc <- src[nzchar(other_taxa_classes)]
    for (i in seq_len(nrow(oc))) {
      cl <- strsplit(strsplit(oc$other_taxa_classes[i], "|", fixed = TRUE)[[1]], ":", fixed = TRUE)
      add(sprintf("Classes left out of %s (rows in the sites):", basename(oc$source[i])), "",
          tab(data.table::data.table(class = vapply(cl, `[`, "", 1), rows = as.numeric(vapply(cl, function(z) z[length(z)], "")))), "")
    }
    if (!all(src$read_completely)) verdict <- c(verdict, "a source file was not read completely")
    if (sum(src$n_wrong_fields)) verdict <- c(verdict, "source lines that could not be read")
  } else add("The reference step has not been run.", "")

  # 3. Cubes against sources
  add("## 3. Site cubes against their sources", "")
  if (!is.null(st)) {
    st[, identical := rows_cube == rows_same & rows_ref == rows_same]
    st[, status := data.table::fcase(
      !in_shapefile, "cube for a site not in the shapefile",
      n_cube_files == 0 & rows_ref == 0, "no cube, no data in the sources",
      n_cube_files == 0, "no cube, but the sources have data",
      identical, "identical to the sources",
      default = "differs from the sources")]
    add(tab(st[, .(sites = .N, rows_in_cubes = sum(rows_cube), rows_in_sources = sum(rows_ref)), by = status][order(-sites)]), "")
    d <- st[status == "differs from the sources"]
    if (nrow(d)) {
      tot <- d[, lapply(.SD, sum), .SDcols = c("missing_before", "missing_from", "extra_inside_site", "extra_outside_site")]
      add(sprintf("Rows in the sources but not in the cube: %s before %d, %s from %d on. Rows in the cube but not in the sources: %s with the cell inside the site, %s with the cell outside it.",
                  fmt(tot$missing_before), from_year, fmt(tot$missing_from), from_year, fmt(tot$extra_inside_site), fmt(tot$extra_outside_site)), "",
          "Sites with the most differing rows (all in `sites.csv`, examples of the rows in `examples.csv`):", "",
          tab(utils::head(d[order(-(missing_before + missing_from + extra_inside_site + extra_outside_site)),
                            .(site_id, rows_cube, rows_ref, missing_before, missing_from, extra_inside_site, extra_outside_site)], 25)), "")
      verdict <- c(verdict, sprintf("cubes that differ from their sources (%s)", fmt(nrow(d))))
    }
    if (any(st$status == "no cube, but the sources have data")) verdict <- c(verdict, "sites with data but no cube")
    if (any(st$status == "cube for a site not in the shapefile")) verdict <- c(verdict, "cubes for sites not in the shapefile")

    add("### Checks on the cubes by themselves", "")
    ic <- st[n_cube_files > 0, .(sites_with_duplicate_rows = sum(dup_cell_year_key > 0), duplicate_rows = sum(dup_cell_year_key),
                                 rows_without_year = sum(missing_year), rows_without_cell = sum(missing_cell),
                                 invalid_cells = sum(invalid_cell), rows_without_species = sum(missing_species))]
    add(tab(ic), "", "Duplicate rows are rows with the same cell, year and species key.", "")
    ks <- st[n_cube_files > 0, .(period = c(paste("before", from_year), paste(from_year, "on")),
                                 old_backbone_keys = c(sum(old_keys_before), sum(old_keys_from)),
                                 new_backbone_keys = c(sum(new_keys_before), sum(new_keys_from)))]
    add("Species keys by period (old backbone: numeric; new backbone: short codes such as 7CYVX; a few new keys",
        "happen to be all digits):", "", tab(ks), "")
    if (ic$duplicate_rows > 0 || ic$rows_without_year > 0 || ic$invalid_cells > 0) verdict <- c(verdict, "problems inside cubes")
  } else add("The compare step has not been run for this tag.", "")

  # 4. Taxonomic and temporal coverage
  add("## 4. Coverage by taxon group and year", "")
  if (!is.null(yg) && nrow(yg)) {
    yg[, period := ifelse(as.numeric(year) < from_year, paste("before", from_year), paste(from_year, "on"))]
    g <- yg[, .(rows_cubes = sum(rows_cube), rows_sources = sum(rows_ref), sites_cubes = data.table::uniqueN(site_id[rows_cube > 0]),
                sites_sources = data.table::uniqueN(site_id[rows_ref > 0])), by = .(group, period)][order(group, period)]
    add("Rows and number of sites with data, per taxon group and period, in the cubes and in the sources:", "", tab(g), "",
        "A group with data before but none from the update period (or the other way round) is not covered",
        "consistently over time, whatever the cubes say.", "")
    y <- yg[as.numeric(year) >= 1990, .(rows_cubes = sum(rows_cube), rows_sources = sum(rows_ref)), by = year][order(year)]
    y[, difference := rows_cubes - rows_sources]
    add("Rows per year from 1990, all sites together:", "", tab(y), "")
  }

  add("## Verdict", "")
  if (length(verdict)) add("Not ready for analysis. To look at:", "", paste0("- ", verdict), "") else
    add("Every cube matches its sources, and no problems were found.", "")
  f <- file.path(out, sprintf("REPORT_%s.md", tag))
  writeLines(L, f, useBytes = TRUE)
  f
}
