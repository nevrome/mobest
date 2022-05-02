#' Summarize the results of the origin search
#'
#' Functions to transform the large origin vector tables to meaningful moving window
# ' summaries
#'
#' @param origin_vectors An object of class \code{mobest_originvectors} as created by
#' \link{determine_origin_vectors}
#' @param origin_summary An object of class \code{mobest_originsummary}
#' as created by \link{summarize_origin_vectors}
#' @param ... (Additional) grouping variables (\code{independent_table_id}, \code{dependent_setting_id},
#' \code{kernel_setting_id}, \code{pred_grid_id}, ...)
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
summarize_origin_vectors <- function(
  origin_vectors, ..., window_start, window_stop, window_width, window_step
) {
  .grouping_var <- rlang::ensyms(...)
  # input check
  checkmate::assert_class(origin_vectors, "mobest_originvectors")
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
          io_run_grouped <- io %>%
            dplyr::group_by(.data[["search_id"]]) %>%
            dplyr::summarise(
              mean_spatial_distance = mean(.data[["ov_dist"]]),
              .groups = "drop"
            )
          if (nrow(io) > 0) {
            io %>%
              dplyr::group_by(
                !!!.grouping_var
              ) %>%
              dplyr::summarise(
                z = mean(c(start, end)),
                undirected_mean_spatial_distance =
                  mean(.data[["ov_dist"]]),
                directed_mean_spatial_distance = sqrt(
                  mean(.data[["field_x"]] - .data[["search_x"]])^2 +
                    mean(.data[["field_y"]] - .data[["search_y"]])^2
                ),
                se_spatial_distance = if (nrow(io_run_grouped) >= 3) {
                  calculate_standard_error(io_run_grouped$mean_spatial_distance)
                } else {
                  Inf
                },
                sd_spatial_distance = if (nrow(io_run_grouped) >= 3) {
                  stats::sd(io_run_grouped$mean_spatial_distance)
                } else {
                  Inf
                }
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
                z = mean(c(start, end)),
                undirected_mean_spatial_distance = NA,
                directed_mean_spatial_distance = NA,
                mean_angle_deg = NA,
                se_spatial_distance = Inf,
                sd_spatial_distance = Inf
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
      usd = tidyr::replace_na(.data[["undirected_mean_spatial_distance"]], 0),
      cumsum_undir_dist = cumsum(.data[["usd"]])
    ) %>%
    dplyr::filter(
      is.na(.data[["undirected_mean_spatial_distance"]])
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
