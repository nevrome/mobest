#' estimate_mobility
#'
#' @param input_grid Either an object of class \code{mobest_input_grid} as
#' produced by \link{search_spatial_origin} or a specific object of class
#' \code{mobest_interpol_grid}. Depending on the input a different mobility
#' estimation algorithm is applied
#' @param mobility_regions test
#'
#' @return test
#'
#' @rdname estimate_mobility
#' @export
estimate_mobility <- function(input_grid, mobility_regions) {
  UseMethod("estimate_mobility")
}

#' @rdname estimate_mobility
#' @export
estimate_mobility.default <- function(input_grid, mobility_regions) {
  stop("x is not an object of class mobest_interpol_grid or mobest_input_grid")
}

#' @rdname estimate_mobility
#' @export
estimate_mobility.mobest_interpol_grid <- function(input_grid, mobility_regions) {

  mob <- input_grid %>% tidyr::pivot_wider(
    id_cols = c("pred_grid_id", "point_id", "x", "y", "z"),
    names_from = "dependent_var_id", values_from = c("mean", "sd")
  ) %>% tidyr::pivot_wider(
    id_cols = c("point_id"),
    names_from = c("pred_grid_id"), values_from = c("mean_C1", "mean_C2", "sd_C1", "sd_C2")
  ) %>%
    dplyr::left_join(
      main_pred_grid
    ) %>%
    dplyr::mutate(
      # delta
      delta_x = delta_x,
      delta_y = delta_y,
      delta_z = delta_z,
      # partial derivatives deriv_x_C*
      deriv_x_C1 = (mean_C1_offset_x - mean_C1_main) / delta_x,
      deriv_x_C2 = (mean_C2_offset_x - mean_C2_main) / delta_x,
      # partial derivatives deriv_y_C*
      deriv_y_C1 = (mean_C1_offset_y - mean_C1_main) / delta_y,
      deriv_y_C2 = (mean_C2_offset_y - mean_C2_main) / delta_y,
      # partial derivatives deriv_z_C*
      deriv_z_C1 = (mean_C1_offset_z - mean_C1_main) / delta_z,
      deriv_z_C2 = (mean_C2_offset_z - mean_C2_main) / delta_z,
      # two directional speeds for each spatial direction x and y
      J_x_C1 = -deriv_z_C1/deriv_x_C1,
      J_y_C1 = -deriv_z_C1/deriv_y_C1,
      J_x_C2 = -deriv_z_C2/deriv_x_C2,
      J_y_C2 = -deriv_z_C2/deriv_y_C2,
      # one combined speed for C1 and C2
      J_x = (J_x_C1 + J_x_C2)/2,
      J_y = (J_y_C1 + J_y_C2)/2,
      # final strength of the speed
      J_final = sqrt(J_x^2 + J_y^2)
    ) %>%
    dplyr::select(
      point_id, x, y, z, sd_C1_main, J_x, J_y, J_final
    ) %>%
    dplyr::mutate(
      angle = unlist(Map(function(x,y) {mobest::vec2deg(c(x,y))}, J_x, J_y)),
      J_final_outlier_removed = ifelse(
        J_final > quantile(J_final, probs = 0.90), NA, J_final
      ),
      J_x_outlier_removed = ifelse(
        abs(J_x) > quantile(abs(J_x), probs = 0.90), NA, J_x
      )
    )

  return(mob)

}

#' @rdname estimate_mobility
#' @export
estimate_mobility.mobest_origin_grid <- function(input_grid, mobility_regions) {

  points_regions <- origin_grid %>%
    dplyr::select(.data[["x"]], .data[["y"]], .data[["point_id"]]) %>%
    unique() %>%
    sf::st_as_sf(
      coords = c("x", "y"),
      crs = sf::st_crs(mobility_regions)
    ) %>%
    sf::st_intersection(mobility_regions) %>%
    sf::st_drop_geometry()

  ori <- origin_grid %>%
    dplyr::left_join(points_regions, by = "point_id")

  speed <- ori %>%
    dplyr::group_by(
      .data[["independent_table_id"]],
      .data[["kernel_setting_id"]],
      .data[["z"]],
      .data[["region_id"]]
    ) %>%
    dplyr::summarise(
      #mean_km_per_decade = mean(.data[["spatial_distance"]])/1000/unique(abs(z-z_origin))*10,
      mean_x_to_origin = mean(.data[["x_to_origin"]]),
      mean_y_to_origin = mean(.data[["y_to_origin"]]),
      mean_km_per_decade = sqrt(mean_x_to_origin^2 + mean_y_to_origin^2)/1000/unique(abs(z-z_origin))*10,
      angle_deg = vec2deg(c(.data[["mean_x_to_origin"]], .data[["mean_y_to_origin"]]))
    ) %>%
    dplyr::ungroup()

  return(speed)

}
