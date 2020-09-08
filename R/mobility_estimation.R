#' estimate_mobility
#'
#' @param interpol_grid_origin test
#' @param mobility_regions test
#'
#' @return test
#'
#' @export
estimate_mobility <- function(interpol_grid_origin, mobility_regions) {

  points_regions <- interpol_grid_origin %>%
    dplyr::select(.data[["x"]], .data[["y"]], .data[["point_id"]]) %>%
    unique() %>%
    sf::st_as_sf(
      coords = c("x", "y"),
      crs = sf::st_crs(mobility_regions)
    ) %>%
    sf::st_intersection(mobility_regions) %>%
    sf::st_drop_geometry()

  ori <- interpol_grid_origin %>%
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
