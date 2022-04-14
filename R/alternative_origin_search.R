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
  # (we don't need to assert properties that are already covered by
  # create_model_grid below)
  checkmate::assert_list(
    independent, types = "mobest_spatiotemporalpositions",
    any.missing = F, min.len = 1, names = "strict"
  )
  checkmate::assert_class(
    dependent, classes = "mobest_observations"
  )
  checkmate::assert_list(
    kernel, types = "mobest_kernelsetting",
    any.missing = F, min.len = 1, names = "strict"
  )
  checkmate::assert_list(
    search_independent, types = "mobest_spatiotemporalpositions",
    any.missing = F, min.len = 1, names = "strict"
  )
  checkmate::assert_class(
    search_dependent, classes = "mobest_observations"
  )
  checkmate::assert_true(all(names(dependent) == names(search_dependent)))
  # prepare data
  search_points <- purrr::map2_dfr(
    names(search_independent), search_independent,
    function(name, x) {
      dplyr::bind_cols(x, search_dependent) %>%
        dplyr::mutate(independent_table_id = name, .before = "id") %>%
        dplyr::mutate(search_z = z - rearview_distance)
    }
  )
  search_field <- search_points$search_z %>% unique() %>%
    purrr::map_dfr(
      function(time_slice) {
        search_space_grid %>%
          geopos_to_spatpos(z = time_slice)
      }
    )
  # construct and run model grid to construct the fields
  model_grid <- create_model_grid(
    independent = independent,
    dependent = dependent,
    kernel = kernel,
    prediction_grid = create_spatpos_multi(
      full_search_field = search_field
    )
  )
  interpol_grid <- run_model_grid(model_grid, quiet = quiet)
  # join search points and fields
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
  # calculate overlap probability
  full_search_table_prob <- full_search_table %>%
    dplyr::mutate(
      probability = dnorm(
        x = measured,
        mean = field_mean,
        sd = field_sd
      )
    )
  # output
  return(full_search_table_prob)
}
