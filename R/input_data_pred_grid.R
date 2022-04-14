#' Create a spatial point grid
#'
#' @param area An object of class \code{sf}. Polygons where the spatial grid should
#' be constructed
#' @param spatial_cell_size Numeric. Size of the output spatial grid cells in the unit of
#' \code{area}. See \code{?sf::st_make_grid} for more info
#' @param temporal_layers Numeric vector. Temporal layers of the requested spatiotemporal grid
#'
#' @return An object of class \code{mobest_spatialpositions}
#'
#' @export
create_prediction_grid <- function(area, spatial_cell_size) {
  # input checks
  checkmate::assert_class(classes = "sf")
  # prepare grid
  space_grid <- area %>%
    sf::st_make_grid(cellsize = spatial_cell_size, what = "centers") %>%
    sf::st_sf() %>%
    sf::st_intersection(area) %>%
    dplyr::mutate(
      x = sf::st_coordinates(.)[,1],
      y = sf::st_coordinates(.)[,2]
    ) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(.data[["x"]], .data[["y"]])
  # compile output
  mobest::create_geopos(
    id = 1:nrow(space_grid),
    x = space_grid$x,
    y = space_grid$y
  )
}
