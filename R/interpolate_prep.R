#' create_prediction_grid
#'
#' @param area test
#' @param spatial_cell_size test
#' @param time_layers test
#'
#' @return test
#'
#' @export
create_prediction_grid <- function(area, spatial_cell_size = 100000, time_layers = seq(-7500, -500, 100)) {

  pred_points_space <- area %>%
    sf::st_make_grid(cellsize = spatial_cell_size, what = "centers") %>%
    sf::st_intersection(area) %>%
    sf::st_coordinates() %>%
    tibble::as_tibble() %>%
    dplyr::rename(x = X, y = Y)

  time_layers <- tibble::tibble(
    z = time_layers
  )

  pred_grid <- pred_points_space %>%
    tidyr::crossing(time_layers) %>%
    dplyr::mutate(
      point_id = 1:nrow(.),
    )

  return(pred_grid)
}
