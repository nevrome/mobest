#' Summarize the results of the origin search
#'
#' Functions to transform the large origin vector tables to meaningful summaries.
#' \code{pack_origin_vectors} yields simple mean vectors,
#' \code{summarize_origin_vectors} a moving window summary through time.
#'
#' @param origin_vectors An object of class \code{mobest_originvectors} as created by
#' \link{determine_origin_vectors} or (for \code{summarize_origin_vectors}) an object
#' of class \code{mobest_originvectorspacked} as created by \link{pack_origin_vectors}
#' @param origin_summary An object of class \code{mobest_originsummary}
#' as created by \link{summarize_origin_vectors}
#' @param ... (Additional) grouping variables (\code{independent_table_id},
#' \code{dependent_setting_id}, \code{kernel_setting_id}, \code{pred_grid_id}, ...)
#' @param window_start Start date of the moving window sequence
#' @param window_stop Stop date of the moving window sequence
#' @param window_width Width of each individual moving window
#' @param window_step Frequency of moving windows. Example:
#' If the first window starts at -3500 and extends until -3200, should the next
#' one start at -3450 or -3400, so with window_step = 50 or window_step = 100?
#'
#' @return \link{summarize_origin_vectors} returns and object of class
#' \code{mobest_originsummary}, \link{find_no_data_windows} then
#' \code{mobest_nodatawindows}.
#'
#' @name origin_summary
NULL

#' @rdname origin_summary
#' @export
pack_origin_vectors <- function(origin_vectors, ...) {
  .grouping_var <- rlang::ensyms(...)
  # input check
  checkmate::assert_class(origin_vectors, "mobest_originvectors")
  # pack vectors
  packed_origin_vectors <- origin_vectors %>%
    dplyr::group_by(
      .data[["search_id"]],
      !!!.grouping_var
    ) %>%
    dplyr::summarise(
      field_x      = mean(.data[["field_x"]]),
      field_y      = mean(.data[["field_y"]]),
      field_z      = mean(.data[["field_z"]]),
      search_x     = mean(.data[["search_x"]]),
      search_y     = mean(.data[["search_y"]]),
      search_z     = mean(.data[["search_z"]]),
      mean_ov_x    = mean(.data[["ov_x"]]),
      mean_ov_y    = mean(.data[["ov_y"]]),
      ov_dist      = sqrt(.data[["mean_ov_x"]]^2 + .data[["mean_ov_y"]]^2),
      ov_dist_se   = calculate_standard_error(sqrt(.data[["ov_x"]]^2 + .data[["ov_y"]]^2)),
      ov_dist_sd   = stats::sd(sqrt(.data[["ov_x"]]^2 + .data[["ov_y"]]^2)),
      ov_angle_deg = vec2deg(c(.data[["mean_ov_x"]], .data[["mean_ov_y"]])),
      .groups = "drop"
    ) %>%
    dplyr::rename(
      ov_x = .data[["mean_ov_x"]],
      ov_y = .data[["mean_ov_y"]]
    )
  # compile output
  packed_origin_vectors %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_originvectorspacked")
}

#' @rdname origin_summary
#' @export
summarize_origin_vectors <- function(
  origin_vectors, ..., window_start, window_stop, window_width, window_step
) {
  .grouping_var <- rlang::ensyms(...)
  # input check
  checkmate::assert(
    checkmate::check_class(origin_vectors, "mobest_originvectorspacked"),
    checkmate::check_class(origin_vectors, "mobest_originvectors")
  )
  checkmate::assert_number(window_start)
  checkmate::assert_number(window_stop, lower = window_start)
  checkmate::assert_number(window_width)
  checkmate::assert_number(window_step)
  # split vector groups
  vector_groups <- origin_vectors %>%
    dplyr::group_split(
      !!!.grouping_var
    )
  # loop through units
  future::plan(future::multisession)
  origin_summary <- furrr::future_map_dfr(
    vector_groups,
    function(vector_group) {
      purrr::map2_df(
        # define moving windows and loop through them
        seq(window_start, window_stop - window_width, window_step),
        seq(window_start + window_width, window_stop, window_step),
        function(start, end) {
          # prepare window data subsets
          io <- dplyr::filter(
            vector_group,
            .data[["search_z"]] >= start,
            .data[["search_z"]] < end
          )
          if (nrow(io) > 0) {
            io %>%
              dplyr::group_by(
                !!!.grouping_var
              ) %>%
              dplyr::summarise(
                z            = mean(c(start, end)),
                ov_dist      = sqrt(mean(.data[["ov_x"]])^2 + mean(.data[["ov_y"]])^2),
                ov_dist_se   = if (dplyr::n() >= 2) {
                                 calculate_standard_error(sqrt(.data[["ov_x"]]^2 + .data[["ov_y"]]^2))
                               } else { Inf },
                ov_dist_sd   = if (dplyr::n() >= 2) {
                                 stats::sd(sqrt(.data[["ov_x"]]^2 + .data[["ov_y"]]^2))
                               } else { Inf },
                ov_angle_deg = vec2deg(c(mean(.data[["ov_x"]]), mean(.data[["ov_x"]])))
              )
            # tibble::tibble(
            #   fraction_smaller_500 = sum(io$spatial_distance < 500) / nrow(io),
            #   fraction_bigger_500 = sum(io$spatial_distance >= 500 & io$spatial_distance < 1000) / nrow(io),
            #   fraction_bigger_1000 = sum(io$spatial_distance >= 1000 & io$spatial_distance < 2000) / nrow(io),
            #   fraction_bigger_2000 = sum(io$spatial_distance >= 2000) / nrow(io)
            # )
          } else {
            vector_group %>%
              dplyr::group_by(
                !!!.grouping_var
              ) %>%
              dplyr::summarise(
                z            = mean(c(start, end)),
                ov_dist      = NA,
                ov_dist_se   = Inf,
                ov_dist_sd   = Inf,
                ov_angle_deg = NA
              )
          }
        }
      )
    }
  )
  origin_summary %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_originsummary")
}

#' @rdname origin_summary
#' @export
find_no_data_windows <- function(origin_summary, ...) {
  .grouping_var <- rlang::ensyms(...)
  # input check
  checkmate::assert_class(origin_summary, "mobest_originsummary")
  # find windows
  window_step <- abs(origin_summary$z[1] - origin_summary$z[2])
  origin_summary %>%
    dplyr::group_by(
      !!!.grouping_var
    ) %>%
    dplyr::mutate(
      usd = tidyr::replace_na(.data[["ov_dist"]], 0),
      cumsum_undir_dist = cumsum(.data[["usd"]])
    ) %>%
    dplyr::filter(
      is.na(.data[["ov_dist"]])
    ) %>%
    dplyr::group_by(
      !!!.grouping_var,
      .data[["cumsum_undir_dist"]]
    ) %>%
    dplyr::summarise(
      min_date_not_covered = min(.data[["z"]]) - window_step,
      max_date_not_covered = max(.data[["z"]]) + window_step,
      .groups = "drop"
    ) %>%
    dplyr::select(-.data[["cumsum_undir_dist"]]) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_nodatawindows")
}

#### helper functions ####

calculate_standard_error <- function(x) { stats::sd(x)/sqrt(length(x)) }
