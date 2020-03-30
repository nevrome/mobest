#' estimate_mobility
#'
#' @param spatial_origin test
#' @param mobility_regions test
#'
#' @return test
#'
#' @export
estimate_mobility <- function(interpol_grid_origin, mobility_regions) {

  points_regions <- interpol_grid_origin %>%
    dplyr::select(x, y, point_id) %>%
    unique() %>%
    sf::st_as_sf(
      coords = c("x", "y"),
      crs = "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"
    ) %>%
    sf::st_intersection(mobility_regions) %>%
    tibble::as_tibble() %>%
    dplyr::select(-geometry)

  ori <- interpol_grid_origin %>%
    dplyr::left_join(points_regions, by = "point_id")

  speed <- ori %>%
    dplyr::group_by(
      independent_table_id, kernel_setting_id, z, region_id
    ) %>%
    dplyr::summarise(
      mean_km_per_decade = mean(spatial_distance)/1000/10,
      mean_x_to_origin = mean(x_to_origin),
      mean_y_to_origin = mean(y_to_origin),
      angle_deg = angle_between_along_360(c(mean_x_to_origin, mean_y_to_origin))
    ) %>%
    dplyr::ungroup()

  return(speed)

}

angle_between_along_360 <- function(x) {
  a <- c(0,1)
  b <- x

  theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) ) * 180 / pi

  if (b[1] < 0) {
    360 - theta
  } else{
    theta
  }
}
