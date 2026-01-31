#' Calculate a Biodiversity Indicator for a Ramsar Site
#'
#' @param cube A data frame or path to a GBIF data cube.
#' @param metric The name of the indicator to calculate (e.g., "obs_richness", "total_occ").
#' @param type The type of indicator ("ts" for time series, "map" for spatial).
#' @param shapefile_path Optional path to the site's WKT/shapefile.
#' @param region The continent or region for the analysis.
#' @param ... Additional arguments passed to b3gbi::compute_indicator_workflow.
#' @return A b3gbi_indicator object.
#'
#' @importFrom b3gbi process_cube compute_indicator_workflow
#' @export
calculate_ramsar_metric <- function(cube,
                                   metric,
                                   type = c("ts", "map"),
                                   shapefile_path = NULL,
                                   region = NULL,
                                   ...) {
  type <- match.arg(type)
  
  # If cube is a path, process it
  if (is.character(cube) && file.exists(cube)) {
    cube <- b3gbi::process_cube(cube, separator = ",")
  }
  
  # Calculate using b3gbi
  indicator_obj <- b3gbi::compute_indicator_workflow(
    data = cube,
    type = metric,
    dim_type = type,
    shapefile_path = shapefile_path,
    region = region,
    ...
  )
  
  return(indicator_obj)
}