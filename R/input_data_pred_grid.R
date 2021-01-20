#' Create a spatiotemporal point grid
#'
#' Creates a prediction grid that can be used in the \link{create_model_grid} +
#' \link{run_model_grid} interpolation workflow
#'
#' @param area An object of class \code{sf}. Polygons where the spatial grid should
#' be constructed
#' @param spatial_cell_size Numeric. Size of the output spatial grid cells in the unit of
#' \code{area}. See \code{?sf::st_make_grid} for more info
#' @param time_layers Numeric vector. Temporal layers of the requested spatiotemporal grid
#'
#' @return Dataframe with columns x, y and z (spatiotemporal coordinates) and a
#' point identifier column id
#'
#' @export
prediction_grid_for_spatiotemporal_area <- function(
  area, regions, spatial_cell_size, temporal_layers) {

  space_grid <- area %>%
    sf::st_make_grid(cellsize = spatial_cell_size, what = "centers") %>%
    sf::st_sf() %>%
    sf::st_intersection(area) %>%
    dplyr::mutate(
      x = sf::st_coordinates(.)[,1],
      y = sf::st_coordinates(.)[,2]
    ) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(
      x, y
    )

  time_grid <- tibble::tibble(
    z = temporal_layers
  )

  pred_grid <- space_grid %>%
    tidyr::crossing(time_grid)

  mobest::create_spatpos(
    id = 1:nrow(pred_grid),
    x = pred_grid$x,
    y = pred_grid$y,
    z = pred_grid$z
  )
}
