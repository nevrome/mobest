#' @rdname search_spatial_origin
#' @export
search_origin <- function(
  independent,
  dependent,
  kernel,
  search_independent,
  search_dependent,
  search_space_grid,
  rearview_distance = 0,
  quiet = F
) {
  # input checks

  # prepare data
  search_points <- purrr::map2_dfr(
    names(search_independent), search_independent,
    function(name, x) {
      dplyr::bind_cols(x, search_dependent) %>%
        dplyr::mutate(independent_table_id = name, .before = "id") %>%
        dplyr::mutate(search_z = z - rearview_distance)
    }
  )
  search_field <- search_points$z %>% unique() %>%
    purrr::map_dfr(
      function(time_slice) {
        search_space_grid %>% dplyr::mutate(z = time_slice)
      }
    )
  model_grid <- mobest::create_model_grid(
    independent = independent,
    dependent = dependent,
    kernel = kernel,
    prediction_grid = create_spatpos_multi(
      full_search_field = search_field
    )
  )
  interpol_grid <- mobest::run_model_grid(model_grid, quiet = T)

  full_search_table <- dplyr::left_join(
    search_points %>%
      tidyr::pivot_longer(
        cols = tidyselect::all_of(names(dependent)),
        names_to = "dependent_var_id",
        values_to = "measured"
      ),
    interpol_grid %>%
      magrittr::set_colnames(paste0("field_", colnames(interpol_grid))),
    by = c(
      "dependent_var_id" = "field_dependent_var_id",
      "search_z" = "field_z"
    )
  )

  hu <- full_search_table %>%
    dplyr::mutate(
      position_probability = dnorm(
        x = measured,
        mean = field_mean,
        sd = field_sd
      )
    )

  gu <- hu %>%
    dplyr::filter(
      independent_table_id == "dating_1",
      id == 3,
      dependent_var_id == "ac1",
      field_independent_table_id == "dating_1",
      field_kernel_setting_id == "kernel_1"
    )

  library(ggplot2)
  gu %>%
    ggplot() +
    geom_raster(
      aes(x = field_x, y = field_y, fill = position_probability)
    )


}
