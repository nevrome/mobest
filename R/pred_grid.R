#' Create a spatiotemporal point grid
#'
#' Creates a prediction grid that can be used in the \link{create_model_grid} +
#' \link{run_model_grid} interpolation workflow#'
#'
#' @param area An object of class \code{sf}. Polygons where the spatial grid should
#' be constructed
#' @param mobility_regions An object of class \code{sf}. Polygons of regions the grid points
#' should be attributed to
#' @param spatial_cell_size Numeric. Size of the spatial grid cells in the unit of
#' \code{area}. See \code{?sf::st_make_grid} for more info
#' @param time_layers Numeric vector. Temporal layers of the requested spatiotemporal
#' grid
#'
#' @return Dataframe with columns x, y and z (spatiotemporal coordinates) and a
#' point identifier column point_id
#'
#' @export
create_prediction_grid <- function(area, mobility_regions, spatial_cell_size = 100000, time_layers = seq(-7500, -500, 100)) {

  point_grid <- area %>%
    sf::st_make_grid(cellsize = spatial_cell_size, what = "centers") %>%
    sf::st_sf() %>%
    sf::st_intersection(area) %>%
    sf::st_join(mobility_regions, join = sf::st_intersects) %>%
    dplyr::mutate(
      x = sf::st_coordinates(.)[,1],
      y = sf::st_coordinates(.)[,2]
    ) %>%
    sf::st_drop_geometry()

  time_grid <- tibble::tibble(
    z = time_layers
  )

  pred_grid <- point_grid %>%
    tidyr::crossing(time_grid) %>%
    dplyr::mutate(
      point_id = 1:nrow(.)
    )

  return(pred_grid)
}
