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
      mean_angle = mean(circular::circular(
        angle_degree, type = "angles", units = "degrees", modulo = "2pi", template = 'geographics'),
        na.rm = T
      )
    ) %>%
    dplyr::ungroup()

  return(speed)

}
