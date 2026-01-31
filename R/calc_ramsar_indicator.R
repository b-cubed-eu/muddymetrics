#' Batch Calculate and Save Ramsar Biodiversity Indicators
#'
#' @description
#' This function iterates through Ramsar sites in country subdirectories,
#' processes their GBIF data cubes, calculates a specific biodiversity
#' indicator using \code{b3gbi}, and saves the resulting plots.
#'
#' @param indicator Character. The indicator to calculate, combined with type
#'   (e.g., "obs_richness_ts", "total_occ_map").
#' @param inputdir Character. Path to the base directory containing country-level
#'   subdirectories with GBIF data cubes.
#' @param maindir Character. Path to the base output directory.
#' @param shapefiledir Character. Path to the base directory containing site WKT files.
#' @param continent Character. The continent name used for b3gbi internal mapping.
#' @param plot_args List. Additional arguments passed to the plot function
#'   (e.g., \code{list(smoothed_trend = FALSE)}).
#' @param ... Additional arguments passed to \code{b3gbi::compute_indicator_workflow}.
#'
#' @return A list containing two nested lists: \code{mean} (mean indicator values
#'   per site) and \code{values} (raw indicator values per site).
#'
#' @importFrom b3gbi process_cube compute_indicator_workflow
#' @importFrom ggplot2 ggsave
#' @importFrom tools file_path_sans_ext
#' @importFrom utils modifyList
#' @export
calc_ramsar_indicator <- function(indicator,
                                  inputdir,
                                  maindir,
                                  shapefiledir,
                                  continent,
                                  plot_args = list(),
                                  ...) {

  # Extract indicator name and type (ts/map)
  if (grepl("_ts$", indicator)) {
    type <- "ts"
    indy <- sub("_ts$", "", indicator)
  } else if (grepl("_map$", indicator)) {
    type <- "map"
    indy <- sub("_map$", "", indicator)
  } else {
    stop("Indicator must end with '_ts' or '_map'")
  }

  countrylist <- list.files(inputdir)
  mean_list <- vector("list", length(countrylist))
  values_list <- vector("list", length(countrylist))
  
  for (i in seq_along(countrylist)) {
    country_name <- countrylist[i]
    country_input_dir <- file.path(inputdir, country_name)
    country_output_dir <- file.path(maindir, country_name)
    
    if (!dir.exists(country_output_dir)) {
      dir.create(country_output_dir, recursive = TRUE)
      message(paste0("Created ", country_name, " output directory: ", country_output_dir))
    }

    sitelist <- list.files(country_input_dir, pattern = "\\.csv$")
    mean_temp <- vector("list", length(sitelist))
    values_temp <- vector("list", length(sitelist))
    
    for (j in seq_along(sitelist)) {
      tryCatch({
        # Robustly extract site name by removing '_data' suffix and extension
        site_file <- sitelist[j]
        sitename <- sub("_data$", "", tools::file_path_sans_ext(site_file))

        shapefilepath <- file.path(shapefiledir, country_name, paste0(sitename, ".wkt"))
        cubepath <- file.path(country_input_dir, site_file)

        cube <- b3gbi::process_cube(cubepath, separator = ",")

        ind_temp <- b3gbi::compute_indicator_workflow(
          data = cube,
          type = indy,
          dim_type = type,
          shapefile_path = shapefilepath,
          ne_scale = "large",
          region = continent,
          include_water = TRUE,
          shapefile_crs = 4326,
          ci_type = "none",
          ...
        )

        values_temp[[j]] <- ind_temp$data$diversity_val
        mean_temp[[j]] <- mean(ind_temp$data$diversity_val, na.rm = TRUE)

        if (type == "ts") {
          # Default TS plot args
          default_plot_args <- list(x = ind_temp, x_expand = 0.1, y_expand = 0.1)
          # Override with user provided args
          final_plot_args <- modifyList(default_plot_args, plot_args)
          plot <- do.call(b3gbi::plot_ts, final_plot_args)
        } else {
          # Default Map plot args
          default_plot_args <- list(x = ind_temp, layers = "lakes")
          # Override with user provided args
          final_plot_args <- modifyList(default_plot_args, plot_args)
          plot <- do.call(b3gbi::plot_map, final_plot_args)
        }

        ggplot2::ggsave(
          filename = file.path(country_output_dir, paste0(sitename, "_", indicator, ".png")),
          plot = plot,
          device = "png",
          height = 4000,
          width = 4000,
          units = "px"
        )

      }, error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("No spatial intersection", msg) ||
            grepl("Error in FUN", msg) ||
            grepl("TopologyException", msg)) {
          message(paste0("Encountered an error for site '", sitename, "': ", msg))
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
