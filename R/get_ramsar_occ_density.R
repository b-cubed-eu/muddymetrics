get_ramsar_occ_density <- function(total_occ_map) {
  x <- total_occ_map
  area <- x$original_bbox |>
    sf::st_area() |>
    units::set_units("km^2")
  occ <- sum(x$data$diversity_val, na.rm = TRUE)
  density <- occ / area
  return(density)
}
