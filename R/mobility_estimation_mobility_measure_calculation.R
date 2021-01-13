#' estimate_mobility
#'
#' @param input_grid Either an object of class \code{mobest_input_grid} as
#' produced by \link{search_spatial_origin} or a specific object of class
#' \code{mobest_interpol_grid}. Depending on the input a different mobility
#' estimation algorithm is applied
#' @param delta_x test
#' @param delta_y test
#' @param delta_z test
#'
#' @return test
#'
#' @rdname estimate_mobility
#' @export
estimate_mobility <- function(input_grid, delta_x, delta_y, delta_z) {
  UseMethod("estimate_mobility")
}

#' @rdname estimate_mobility
#' @export
estimate_mobility.default <- function(input_grid, delta_x, delta_y, delta_z) {
  stop("x is not an object of class mobest_interpol_grid or mobest_input_grid")
}

#' @rdname estimate_mobility
#' @export
estimate_mobility.mobest_origin_grid <- function(input_grid, delta_x, delta_y, delta_z) {

  speed <- input_grid %>%
    dplyr::mutate(
      speed_km_per_decade = .data[["spatial_distance"]]/1000/unique(abs(.data[["z"]]-.data[["z_origin"]]))*10
    )

  ## mean by region
  # speed <- ori %>%
  #   dplyr::group_by(
  #     .data[["independent_table_id"]],
  #     .data[["kernel_setting_id"]],
  #     .data[["z"]],
  #     .data[["region_id"]]
  #   ) %>%
  #   dplyr::summarise(
  #     #mean_km_per_decade = mean(.data[["spatial_distance"]])/1000/unique(abs(z-z_origin))*10,
  #     mean_x_to_origin = mean(.data[["x_to_origin"]]),
  #     mean_y_to_origin = mean(.data[["y_to_origin"]]),
  #     mean_km_per_decade = sqrt(mean_x_to_origin^2 + mean_y_to_origin^2)/1000/unique(abs(z-z_origin))*10,
  #     angle_deg = vec2deg(c(.data[["mean_x_to_origin"]], .data[["mean_y_to_origin"]]))
  #   ) %>%
  #   dplyr::ungroup()

  return(speed)

}
