#' Generate GBIF Predicates from WKT and Filters
#'
#' @param wkt A WKT string defining the geometry.
#' @param filters A list of filters to apply (e.g., hasGeospatialIssue, basisOfRecord).
#' @return A GBIF predicate object.
#' @export
get_gbif_predicates <- function(wkt, filters) {
  # Start with the geometry predicate
  predicates <- list(rgbif::pred("geometry", wkt))
  
  # Add other filters dynamically
  for (n in names(filters)) {
    val <- filters[[n]]
    if (length(val) == 1) {
      predicates[[length(predicates) + 1]] <- rgbif::pred(n, val)
    } else {
      predicates[[length(predicates) + 1]] <- rgbif::pred_in(n, val)
    }
  }
  
  # Combine with AND
  do.call(rgbif::pred_and, predicates)
}
