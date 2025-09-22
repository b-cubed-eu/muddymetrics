calc_ramsar_indicator <- function(indicator,
                                  inputdir,
                                  maindir,
                                  shapefiledir,
                                  continent) {

  type = substring(indicator, nchar(indicator) - 1, nchar(indicator))
  type <- ifelse(type == "ap", "map", "ts")
  indy <- substring(
    indicator, 1, ifelse(
      type == "ts", nchar(indicator) - 3, nchar(indicator) - 4
    )
  )

  countrylist <- list.files(inputdir)
  mean_list <- vector("list", length(countrylist))
  values_list <- vector("list", length(countrylist))
  for (i in 1:length(countrylist)) {

    countrydir <- paste0(maindir, "/", countrylist[i])
    if (!dir.exists(countrydir)) {
      dir.create(countrydir, recursive = TRUE)
      message(paste0("Created ", countrylist[i], " output directory:",
                     countrydir))
    }

    sitelist <- list.files(paste0(inputdir, "/", countrylist[i]))

    mean_temp <- vector("list", length(sitelist))
    values_temp <- vector("list", length(sitelist))
    for (j in 1:length(sitelist)) {

      tryCatch({

        sitename <- substring(sitelist[j], 1, nchar(sitelist[j]) - 9)

        shapefilepath <- paste0(
          shapefiledir, "/", countrylist[i], "/", sitename, ".wkt"
        )

        cube <- process_cube(paste0(
          inputdir, "/", countrylist[i], "/", sitelist[j]), separator = ","
        )

        ind_temp <- compute_indicator_workflow(
          data = cube,
          type = indy,
          dim_type = type,
          shapefile_path = shapefilepath,
          ne_scale = "large",
          region = continent,
          include_water = TRUE,
          shapefile_crs = 4326,
          ci_type = "none"
        )

        values_temp[[j]] <- ind_temp$data$diversity_val
        mean_temp[[j]] <- mean(ind_temp$data$diversity_val, na.rm = TRUE)

        if (type == "ts") {
          plot <- plot(ind_temp, x_expand = 0.1, y_expand = 0.1)
        } else {
          plot <- plot(ind_temp, layers = "lakes")
        }

        ggplot2::ggsave(filename = paste0(
          maindir, "/", countrylist[i], "/", sitename, "_", indicator, ".png"
          ), plot, device = "png", height = 4000, width = 4000, units = "px")

      }, error = function(e) {
        if (grepl("No spatial intersection between map data and grid.", e) ||
            grepl("Error in FUN", e) ||
            grepl("Loop", e) ||
            grepl("TopologyException", e)) {
          message(
            paste0(
              "Encountered a geometry error during intersection. This may be ",
              "due to coordinate mismatches."
            )
          )
        } else {
          stop(e)
        }
      })

    }
    names(mean_temp) <- sitelist
    names(values_temp) <- sitelist
    mean_list[[i]] <- mean_temp
    values_list[[i]] <- values_temp
  }

  names(values_list) <- countrylist
  names(mean_list) <- countrylist

  return(list("mean" = mean_list, "values" = values_list))

}
