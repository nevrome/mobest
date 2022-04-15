#' @rdname search_spatial_origin
#' @export
locate <- function(
  independent,
  dependent,
  kernel,
  search_independent,
  search_dependent,
  search_space_grid,
  search_time = 0,
  search_time_mode = "relative",
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
  checkmate::assert_numeric(
    search_time,
    finite = TRUE, any.missing = FALSE, min.len = 1, unique = TRUE
  )
  checkmate::assert_choice(
    search_time_mode, choices = c("relative", "absolute")
  )
  checkmate::assert_true(all(names(dependent) == names(search_dependent)))
  # prepare data
  search_points <- purrr::map2_dfr(
    names(search_independent), search_independent,
    function(name, x) {
      dplyr::bind_cols(x, search_dependent) %>%
        dplyr::mutate(independent_table_id = name, .before = "id") %>%
        tidyr::crossing(tibble::tibble(search_time = search_time)) %>%
        dplyr::mutate(
          search_z =
            if (search_time_mode == "relative") {
              .data[["z"]] + .data[["search_time"]]
            } else if (search_time_mode == "absolute") {
              .data[["search_time"]]
            }
        ) %>%
        dplyr::select(-.data[["search_time"]])
    }
  )
  search_fields <- purrr::map(
    search_points$search_z %>% unique(),
    function(time_slice) {
      search_space_grid %>% geopos_to_spatpos(z = time_slice)
    }
  ) %>% magrittr::set_names(., paste("time_slice", 1:length(.), sep = "_"))
  # construct and run model grid to construct the fields
  model_grid <- create_model_grid(
    independent = independent,
    dependent = dependent,
    kernel = kernel,
    prediction_grid = do.call(create_spatpos_multi, search_fields)
  )
  if (!quiet) { message("Constructing search fields") }
  interpol_grid <- run_model_grid(model_grid, quiet = quiet)
  # join search points and fields
  if (!quiet) { message("Compiling full search table") }
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
  if (!quiet) { message("Calculating probabilities") }
  full_search_table_prob <- full_search_table %>%
    dplyr::mutate(
      probability = dnorm(
        x = measured,
        mean = field_mean,
        sd = field_sd
      )
    )
  # output
  full_search_table_prob %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_locateoverview") %>%
    return()
}

#' @rdname search_spatial_origin
#' @export
multiply_dependent_probabilities <- function(locate_overview, omit_dependent_details = T) {
  # input checks
  checkmate::assert_class(locate_overview, classes = "mobest_locateoverview")
  # data transformation
  locate_summary <- locate_overview %>%
    tidyr::pivot_wider(
      names_from = "dependent_var_id",
      values_from = c(
        "measured",
        "field_dsx", "field_dsy", "field_dt", "field_g",
        "field_mean", "field_sd",
        "probability",
      )
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      probability_product = purrr::reduce(
        dplyr::c_across(tidyselect::starts_with("probability")),
        function(x, y) { x*y }
      )
    ) %>%
    dplyr::ungroup()
  # prepare output
  if (omit_dependent_details) {
    locate_summary %>%
      dplyr::select(-tidyselect::ends_with(
        locate_overview$dependent_var_id %>% unique
      ))
  } else {
    locate_summary
  }
}
