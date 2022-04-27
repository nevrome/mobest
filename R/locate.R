#' @rdname search_spatial_origin
#' @export
locate <- function(
  independent,
  dependent,
  kernel,
  search_independent,
  search_dependent,
  search_dependent_error = NULL,
  search_space_grid,
  search_time = 0,
  search_time_mode = "relative",
  quiet = F
) {
  locate_multi(
    independent = create_spatpos_multi(i = independent),
    dependent = dependent,
    kernel = create_kernset_multi(k = kernel),
    search_independent = create_spatpos_multi(si = search_independent),
    search_dependent = search_dependent,
    search_dependent_error = search_dependent_error,
    search_space_grid = search_space_grid,
    search_time = search_time,
    search_time_mode = search_time_mode,
    quiet = F
  ) %>%
    dplyr::select(
      -.data[["independent_table_id"]],
      -.data[["field_independent_table_id"]],
      -.data[["field_kernel_setting_id"]]
    )
}

#' @rdname search_spatial_origin
#' @export
locate_multi <- function(
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
  checkmate::assert_class(search_independent, "mobest_spatiotemporalpositions_multi")
  checkmate::assert(
    checkmate::check_class(
      search_dependent, classes = "mobest_observations_multi"
    ),
    checkmate::check_class(
      search_dependent, classes = "mobest_observations_with_error_multi"
    )
  )
  checkmate::assert_class(search_space_grid, "mobest_spatialpositions")
  checkmate::assert_numeric(
    search_time,
    finite = TRUE, any.missing = FALSE, min.len = 1, unique = TRUE
  )
  checkmate::assert_choice(
    search_time_mode, choices = c("relative", "absolute")
  )
  check_compatible_multi(search_independent, search_dependent, check_df_nrow_equal)
  check_compatible_multi(dependent, search_dependent, check_names_equal, ignore_sd_cols = T)
  checkmate::assert_true(setequal(names(independent), names(search_independent)))
  checkmate::assert_true(setequal(names(dependent), names(search_dependent)))
  # prepare data
  search_points <- tidyr::crossing(
    search_independent = search_independent %>%
      purrr::map2(names(.), ., function(n, x) { x$independent_table_id <- n; x}),
    search_dependent = search_dependent %>%
      purrr::map2(names(.), ., function(n, x) { x$dependent_setting_id <- n; x})
  ) %>% tidyr::unnest(
    cols = c("search_independent", "search_dependent")
  ) %>%
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
        cols = -c(
          "id", "x", "y", "z", "independent_table_id", "dependent_setting_id", "search_z"
        ),
        names_to = "dependent_var_id",
        values_to = "intermediate_value"
      ) %>%
      tidyr::separate(
        col = "dependent_var_id",
        into = c("dependent_var_id", "dep_var_type"),
        sep = "_(?=sd$)", # only split, if the string ends with "_sd"
        extra = "merge",
        fill = "right"
      ) %>%
      dplyr::mutate(dep_var_type = tidyr::replace_na(dep_var_type, "measured")) %>%
      tidyr::pivot_wider(
        names_from = "dep_var_type",
        values_from = "intermediate_value"
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
  #return(full_search_table)
  if (is.null(search_dependent_error)) {
    full_search_table_prob <- full_search_table %>%
      dplyr::mutate(
        probability = dnorm(
          x = .data[["measured"]],
          mean = .data[["field_mean"]],
          sd = .data[["field_sd"]]
        )
      ) %>%
      dplyr::mutate(
        sd = NA_real_,
        .after = "measured"
      )
  } else {
    int_f <- function(x, mu1, mu2, sd1, sd2) {
      f1 <- dnorm(x, mean = mu1, sd = sd1)
      f2 <- dnorm(x, mean = mu2, sd = sd2)
      pmin(f1, f2)
    }
    full_search_table_prob <- full_search_table %>%
      dplyr::mutate(
        .,
        probability = purrr::pmap_dbl(
          list(.data[["measured"]], .data[["field_mean"]], .data[["sd"]], .data[["field_sd"]]),
          function(mu1, mu2, sd1, sd2) {
            res <- try(
              integrate(
                int_f, -Inf, Inf,
                mu1 = mu1, mu2 = mu2, sd1 = sd1, sd2 = sd2#,
                #rel.tol = 1e-10
              ),
              silent = TRUE
            )
            if(inherits(res ,'try-error')){
              #warning(as.vector(res))
              return(0)
            } else {
              return(res$value)
            }
          }
        )
      )
  }
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
        "measured", "sd",
        "field_dsx", "field_dsy", "field_dt", "field_g",
        "field_mean", "field_sd",
        "probability",
      )
    ) %>%
    dplyr::mutate(
      probability_product = dplyr::select(., tidyselect::starts_with("probability")) %>%
        as.matrix() %>%
        apply(1, prod)
    )
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
