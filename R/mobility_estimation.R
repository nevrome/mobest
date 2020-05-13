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
      crs = "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"
    ) %>%
    sf::st_intersection(mobility_regions) %>%
    tibble::as_tibble() %>%
    dplyr::select(-.data[["geometry"]])

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
      mean_km_per_decade = mean(.data[["spatial_distance"]])/1000/10,
      mean_x_to_origin = mean(.data[["x_to_origin"]]),
      mean_y_to_origin = mean(.data[["y_to_origin"]]),
      angle_deg = vec2deg(c(.data[["mean_x_to_origin"]], .data[["mean_y_to_origin"]]))
    ) %>%
    dplyr::ungroup()

  return(speed)

}
