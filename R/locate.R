#' Calculate a similarity/"origin" probability field for samples of interest
#'
#' Create a genetic ancestry field and determine similarity probabilities with it.
#' \link{locate_multi} handles permutations of input parameters, \link{locate} is
#' a simplified interface.
#'
#' @param independent An object of class \code{mobest_spatiotemporalpositions_multi}
#' or \code{mobest_spatiotemporalpositions_multi} as created by \link{create_spatpos}
#' or \link{create_spatpos_multi}. Spatiotemporal input point positions to inform the
#' ancestry field.
#' @param dependent An object of class \code{mobest_observations} or
#' \code{mobest_observations_multi} as created by \link{create_obs} or
#' \link{create_obs_multi}. Dependent variables that should be interpolated for the
#' ancestry field.
#' @param kernel An object of class \code{mobest_kernelsetting} or
#' \code{mobest_kernelsetting_multi} as created by \link{create_kernset} or
#' \link{create_kernset_multi}. Kernel parameter settings for the ancestry fields.
#' @param search_independent An object of class \code{mobest_spatiotemporalpositions_multi}
#' or \code{mobest_spatiotemporalpositions_multi} as created by \link{create_spatpos}
#' or \link{create_spatpos_multi}. Spatiotemporal input point positions of the samples
#' of interest for which the similarity probability should be calculated. Must have
#' the same names as \code{independent}.
#' @param search_dependent An object of class \code{mobest_observations} or
#' \code{mobest_observations_multi} as created by \link{create_obs} or
#' \link{create_obs_multi}. Dependent variables of the samples of interest for which
#' the similarity probability should be calculated. Must have the same names as
#' \code{dependent}.
#' @param search_space_grid An object of class \code{mobest_spatialpositions}
#' as for example created by \link{create_prediction_grid}. Prediction grid positions
#' for which the similarity probabilities should be printed.
#' @param search_time Numeric vector. Time slices of the predition grid.
#' @param search_time_mode Character choice. One of "relative" or "absolute". Should
#' the search time be relative to the dating of the samples of interest or just absolute
#' points in time? The default for \code{search_time} and \code{search_time_mode} is
#' "relative" and 0, which causes the probabilities to be calculated for the exact dating
#' of the samples of interest.
#' @param normalize Logical. Should the output probability distribution for one permutation
#' be normalized?
#' @param quiet Logical. Should a progress indication be printed?
#' @param locate_overview An object of class \code{mobest_locateoverview} as created by
#' \link{locate} and \link{locate_multi}.
#'
#' @name locate
NULL

#' @rdname locate
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
  normalize = T,
  quiet = F
) {
  locate_multi(
    independent = create_spatpos_multi(i = independent),
    dependent = create_obs_multi(d = dependent),
    kernel = create_kernset_multi(k = kernel),
    search_independent = create_spatpos_multi(i = search_independent),
    search_dependent = create_obs_multi(d = search_dependent),
    search_space_grid = search_space_grid,
    search_time = search_time,
    search_time_mode = search_time_mode,
    normalize = normalize,
    quiet = quiet
  )# %>%
    # dplyr::select(
    #   -.data[["independent_table_id"]],
    #   -.data[["dependent_setting_id"]],
    #   -.data[["kernel_setting_id"]],
    #   -.data[["pred_grid_id"]]
    # )
}

#' @rdname locate
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
  normalize = T,
  quiet = F
) {
  # input checks
  # (we don't need to assert properties that are already covered by
  # create_model_grid below)
  checkmate::assert_class(search_independent, "mobest_spatiotemporalpositions_multi")
  checkmate::assert_class(search_dependent, classes = "mobest_observations_multi")
  checkmate::assert_class(search_space_grid, "mobest_spatialpositions")
  checkmate::assert_numeric(search_time, finite = TRUE, any.missing = FALSE, min.len = 1, unique = TRUE)
  checkmate::assert_choice(search_time_mode, choices = c("relative", "absolute"))
  checkmate::assert_true(setequal(names(independent), names(search_independent)))
  checkmate::assert_true(setequal(names(dependent), names(search_dependent)))
  check_compatible_multi(search_independent, search_dependent, check_df_nrow_equal)
  check_compatible_multi(dependent, search_dependent, check_names_equal, ignore_sd_cols = T)
  # construct search point permutations
  search_points <- tidyr::crossing(
    search_independent = search_independent %>%
      purrr::map2(names(.), ., function(n, x) {x$independent_table_id <- n; x}),
    search_dependent = search_dependent %>%
      purrr::map2(names(.), ., function(n, x) {x$dependent_setting_id <- n; x})
  ) %>% tidyr::unnest(
    cols = c("search_independent", "search_dependent")
  ) %>%
    dplyr::rename_with(
      function(x) { paste0("search_", x) },
      tidyselect::any_of(c("id", "x", "y", "z", "independent_table_id"))
    ) %>%
    tidyr::crossing(tibble::tibble(search_time = search_time)) %>%
    dplyr::mutate(field_z =
      if (search_time_mode == "relative") {
        .data[["search_z"]] + .data[["search_time"]]
      } else if (search_time_mode == "absolute") {
        .data[["search_time"]]
      }
    ) %>%
    tidyr::pivot_longer(
      cols = -c(
        "search_id", "search_x", "search_y", "search_z",
        "search_independent_table_id",
        "dependent_setting_id", "field_z", "search_time"
      ),
      names_to = "dependent_var_id",
      values_to = "search_measured"
    )
  # construct search fields (point rasters for the search)
  unique_time_slices <- tibble::tibble(
    field_z = search_points$field_z %>% unique(),
    pred_grid_id = paste0("time_slice_", seq_along(.data[["field_z"]]))
  )
  search_fields <- purrr::map(
    unique_time_slices$field_z,
    function(time_slice) {search_space_grid %>% geopos_to_spatpos(z = time_slice)}
  ) %>% magrittr::set_names(., unique_time_slices$pred_grid_id)
  # construct model grid for the fields with all (!) permutations
  model_grid <- create_model_grid(
    independent = independent,
    dependent = dependent,
    kernel = kernel,
    prediction_grid = do.call(create_spatpos_multi, search_fields)
  ) %>%
    dplyr::left_join(unique_time_slices, by = "pred_grid_id")
  # remove unnecessary permutations
  model_grid_filtered <- dplyr::semi_join(
    model_grid, search_points,
    by = c(
      "independent_table_id" = "search_independent_table_id",
      "dependent_setting_id",
      "dependent_var_id",
      "field_z"
    )
  ) %>% dplyr::select(-.data[["field_z"]])
  # run model grid to create search fields
  if (!quiet) { message("Constructing search fields") }
  interpol_grid <- run_model_grid(model_grid_filtered, quiet = quiet) %>%
    dplyr::rename_with(
      function(x) { paste0("field_", x) },
      tidyselect::any_of(c("id", "geo_id", "x", "y", "z", "mean", "sd"))
    )
  # join search fields and search points
  if (!quiet) { message("Compiling full search table") }
  full_search_table <- dplyr::left_join(
    interpol_grid, search_points,
    by = c(
      "independent_table_id" = "search_independent_table_id",
      "dependent_setting_id",
      "dependent_var_id",
      "field_z"
    )
  )
  # calculate overlap probabilities
  if (!quiet) { message("Calculating probabilities") }
  full_search_table_prob <- full_search_table %>%
    dplyr::mutate(
      probability = stats::dnorm(
        x = .data[["search_measured"]],
        mean = .data[["field_mean"]],
        sd = .data[["field_sd"]]
      )
    )
  # normalise output probabilities per permutation
  if (normalize) {
    full_search_table_prob <- full_search_table_prob %>%
      dplyr::group_by(
        .data[["independent_table_id"]],
        .data[["dependent_setting_id"]],
        .data[["dependent_var_id"]],
        .data[["kernel_setting_id"]],
        .data[["pred_grid_id"]],
        .data[["search_id"]],
        .data[["search_z"]],
        .data[["search_time"]]
      ) %>%
      dplyr::mutate(
        probability = .data[["probability"]]/sum(.data[["probability"]])
      ) %>%
      dplyr::ungroup()
  }
  # output
  full_search_table_prob %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_locateoverview") %>%
    return()
}

#' @param omit_dependent_details Logical. If TRUE, removes additional, dependent-wise
#' columns.
#'
#' @rdname locate
#' @export
multiply_dependent_probabilities <- function(locate_overview, normalize = T, omit_dependent_details = T) {
  # input checks
  checkmate::assert_class(locate_overview, classes = "mobest_locateoverview")
  # data transformation
  locate_summary <- locate_overview %>%
    tidyr::pivot_wider(
      names_from = "dependent_var_id",
      values_from = tidyselect::any_of(c(
        "search_measured", "search_sd",
        "dsx", "dsy", "dt", "g",
        "field_mean", "field_sd",
        "probability"
      ))
    ) %>%
    dplyr::mutate(
      probability = dplyr::select(., tidyselect::starts_with("probability")) %>%
        as.matrix() %>%
        apply(1, prod)
    )
  # normalise output probabilities per permutation
  if (normalize) {
    locate_summary <- locate_summary %>%
      dplyr::group_by(
        .data[["independent_table_id"]],
        .data[["dependent_setting_id"]],
        .data[["kernel_setting_id"]],
        .data[["pred_grid_id"]],
        .data[["search_id"]],
        .data[["search_z"]],
        .data[["search_time"]]
      ) %>%
      dplyr::mutate(
        probability = .data[["probability"]]/sum(.data[["probability"]])
      ) %>%
      dplyr::ungroup()
  }
  # prepare output
  if (omit_dependent_details) {
    losum <- locate_summary %>%
      dplyr::select(-tidyselect::ends_with(
        locate_overview$dependent_var_id %>% unique
      ))
  } else {
    losum <- locate_summary
  }
  losum %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_locateproduct")
}

#' @param locate_product An object of class \code{mobest_locateproduct}.
#' @param ... (Additional) grouping variables (\code{independent_table_id},
#' \code{dependent_setting_id}, \code{kernel_setting_id}, \code{pred_grid_id}, ...)
#' @param folding_operation Function. Folding operation that should be applied to
#' the probabilities in a group. Default: sum
#'
#' @rdname locate
#' @export
fold_probabilities_per_group <- function(
  locate_product,
  ...,
  folding_operation = function(x) { sum(x, na.rm = T) },
  normalize = T
) {
  .grouping_var <- rlang::ensyms(...)
  # input checks
  checkmate::assert_class(locate_product, classes = "mobest_locateproduct")
  # data transformation
  locate_summary <- locate_product %>%
    dplyr::group_by(
      !!!.grouping_var,
      .data[["search_id"]],
      .data[["field_id"]],
      .data[["search_z"]]
    ) %>%
    dplyr::summarise(
      search_x = dplyr::first(.data[["search_x"]]),
      search_y = dplyr::first(.data[["search_y"]]),
      field_x = dplyr::first(.data[["field_x"]]),
      field_y = dplyr::first(.data[["field_y"]]),
      field_z = dplyr::first(.data[["field_z"]]),
      probability = folding_operation(.data[["probability"]]),
      .groups = "drop"
    )
  # normalise output probabilities per permutation
  if (normalize) {
    locate_summary <- locate_summary %>%
      dplyr::group_by(
        !!!.grouping_var,
        .data[["search_id"]],
        .data[["search_z"]]
      ) %>%
      dplyr::mutate(
        probability = .data[["probability"]]/sum(.data[["probability"]])
      ) %>%
      dplyr::ungroup()
  }
  # prepare output
  locate_summary %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_locatefold")
}
