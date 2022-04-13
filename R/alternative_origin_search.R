#' @rdname search_spatial_origin
#' @export
search_origin <- function(
  independent,
  dependent,
  kernel,
  search_independent,
  search_dependent,
  search_area,
  search_resolution,
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

  search_points$z %>% unique() %>% length()
  prediction_grid_for_spatiotemporal_area(

  )

  model_grid <- mobest::create_model_grid(
    independent = uncertain_positions,
    dependent = observations,
    kernel = mobest::create_kernset_multi(
      kernel_1 = mobest::create_kernset(
        ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
        ac2 = mobest::create_kernel(1000000, 1000000, 200, 0.1)
      ),
      kernel_2 = mobest::create_kernset(
        ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
        ac2 = mobest::create_kernel(1000000, 1000000, 250, 0.1)
      )
    ),
    prediction_grid = mobest::create_spatpos_multi(
      pred_grid_1 = expand.grid(
        x = seq(100000, 1000000, 100000),
        y = seq(100000, 1000000, 100000),
        z = seq(-5500, -3000, 500)
      ) %>% { mobest::create_spatpos(id = 1:nrow(.), x = .$x, y = .$y, z = .$z) }
    )
  )


  search_points <- search_independent

}
