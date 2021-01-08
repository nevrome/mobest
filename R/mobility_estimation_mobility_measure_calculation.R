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
estimate_mobility.mobest_interpol_grid <- function(input_grid, delta_x, delta_y, delta_z) {

  mob <- input_grid %>% tidyr::pivot_wider(
    id_cols = c("pred_grid_id", "point_id", "x", "y", "z", "kernel_setting_id", "independent_table_id"),
    names_from = "dependent_var_id", values_from = c("mean", "sd")
  ) %>% tidyr::pivot_wider(
    id_cols = c("point_id", "kernel_setting_id", "independent_table_id"),
    names_from = c("pred_grid_id"), values_from = c("mean_C1", "mean_C2", "sd_C1", "sd_C2")
  ) %>% dplyr::left_join(
    input_grid %>%
      dplyr::filter(.data[["pred_grid_id"]] == "main") %>%
      dplyr::select(c("point_id", "x", "y", "z", "region_id")) %>%
      unique(),
    by = "point_id"
  ) %>%
  dplyr::mutate(
    # delta
    delta_x = delta_x,
    delta_y = delta_y,
    delta_z = delta_z,
    # partial derivatives deriv_x_C*
    deriv_x_C1 = (.data[["mean_C1_offset_x"]] - .data[["mean_C1_main"]]) / .data[["delta_x"]],
    deriv_x_C2 = (.data[["mean_C2_offset_x"]] - .data[["mean_C2_main"]]) / .data[["delta_x"]],
    # partial derivatives deriv_y_C*
    deriv_y_C1 = (.data[["mean_C1_offset_y"]] - .data[["mean_C1_main"]]) / .data[["delta_y"]],
    deriv_y_C2 = (.data[["mean_C2_offset_y"]] - .data[["mean_C2_main"]]) / .data[["delta_y"]],
    # partial derivatives deriv_z_C*
    deriv_z_C1 = (.data[["mean_C1_offset_z"]] - .data[["mean_C1_main"]]) / .data[["delta_z"]],
    deriv_z_C2 = (.data[["mean_C2_offset_z"]] - .data[["mean_C2_main"]]) / .data[["delta_z"]],
    # two directional speeds for each spatial direction x and y
    speed_x_C1 = -.data[["deriv_z_C1"]]/.data[["deriv_x_C1"]],
    speed_y_C1 = -.data[["deriv_z_C1"]]/.data[["deriv_y_C1"]],
    speed_x_C2 = -.data[["deriv_z_C2"]]/.data[["deriv_x_C2"]],
    speed_y_C2 = -.data[["deriv_z_C2"]]/.data[["deriv_y_C2"]],
    # one combined speed for C1 and C2
    speed_x = (.data[["speed_x_C1"]] + .data[["speed_x_C2"]])/2,
    speed_y = (.data[["speed_y_C1"]] + .data[["speed_y_C2"]])/2,
    # final strength of the speed
    speed_final = sqrt(.data[["speed_x"]]^2 + .data[["speed_y"]]^2),
    # angle
    angle = unlist(Map(function(x,y) {mobest::vec2deg(c(x,y))}, .data[["speed_x"]], .data[["speed_y"]])),
    # rescaling to km/decade
    # speed_x = (speed_x/1000)*10,
    # speed_y = (speed_y/1000)*10,
    speed_km_per_decade = (speed_final/1000)*10
  ) %>%
  dplyr::select(
    "kernel_setting_id", "independent_table_id", "point_id",
    "x", "y", "z", "region_id", "sd_C1_main", "speed_km_per_decade", "angle"
  )

  return(mob)
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
