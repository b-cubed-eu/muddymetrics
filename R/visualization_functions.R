#' Plot and Save a Ramsar Biodiversity Indicator
#'
#' @param indicator_obj A b3gbi_indicator object.
#' @param filename The path to save the plot (including extension).
#' @param device The device to use (default: "png").
#' @param width Width of the plot (default: 4000).
#' @param height Height of the plot (default: 4000).
#' @param units Units for width and height (default: "px").
#' @param ... Additional arguments passed to b3gbi::plot.
#' @return The ggplot object (invisibly).
#' @export
save_ramsar_plot <- function(indicator_obj,
                            filename,
                            device = "png",
                            width = 4000,
                            height = 4000,
                            units = "px",
                            ...) {
  
  # Determine plot type based on object class
  if (inherits(indicator_obj, "indicator_ts")) {
    p <- b3gbi::plot_ts(indicator_obj, ...)
  } else if (inherits(indicator_obj, "indicator_map")) {
    p <- b3gbi::plot_map(indicator_obj, ...)
  } else {
    # Generic plot if class is unknown but b3gbi handles it
    p <- plot(indicator_obj, ...)
  }
  
  # Save using ggplot2
  ggplot2::ggsave(
    filename = filename,
    plot = p,
    device = device,
    width = width,
    height = height,
    units = units
  )
  
  return(invisible(p))
}
