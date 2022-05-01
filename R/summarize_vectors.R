#' Summarize the results of the origin search
#'
#' Functions to transform the large origin grid objects to meaningful summary datasets
#'
#' @param origin_grid An object of class \code{mobest_origingrid} as created by
#' \link{search_spatial_origin}
#' @param moving_origin_grid An object of class \code{mobest_movingorigingrid}
#' as created by \link{average_origin_moving_window}
#' @param window_start Start date of the moving window sequence
#' @param window_stop Stop date of the moving window sequence
#' @param window_width Width of each individual moving window
#' @param window_step Frequency of moving windows. Example:
#' If the first window starts at -3500 and extends until -3200, should the next
#' one start at -3450 or -3400, so with window_step = 50 or window_step = 100?
#'
#' @return Different data products: \code{mobest_meanorigingrid},
#' \code{mobest_movingorigingrid} or \code{mobest_origingridnodatawindows}
#'
#' @name average_origin
NULL

#' @rdname average_origin
#' @export
summarize_origin_vectors <- function(
  origin_vectors, ..., window_start, window_stop, window_width, window_step
) {
  .grouping_var <- rlang::ensyms(...)
  # input check
  checkmate::assert_class(origin_vectors, "mobest_originvectors")
  # split vector groups
  vector_groups <- origin_vectors %>%
    dplyr::group_split(
      !!!.grouping_var
    )
  # loop through units
  future::plan(future::multisession)
  vector_summary <- furrr::future_map_dfr(
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
                  se(io_run_grouped$mean_spatial_distance)
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
                region_id = region,
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
  vector_summary %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_originsummary")
}

se <- function(x) stats::sd(x)/sqrt(length(x))

#' @rdname average_origin
#' @export
no_data_windows <- function(moving_origin_grid, window_step) {
  moving_origin_grid %>%
    dplyr::group_by(.data[["region_id"]]) %>%
    dplyr::mutate(
      usd = tidyr::replace_na(.data[["undirected_mean_spatial_distance"]], 0),
      cumsum_undir_dist = cumsum(.data[["usd"]])
    ) %>%
    dplyr::filter(
      is.na(.data[["undirected_mean_spatial_distance"]])
    ) %>%
    dplyr::group_by(.data[["region_id"]], .data[["cumsum_undir_dist"]]) %>%
    dplyr::summarise(
      min_date_not_covered = min(.data[["z"]]) - window_step,
      max_date_not_covered = max(.data[["z"]]) + window_step,
      .groups = "drop"
    ) %>%
    dplyr::select(-.data[["cumsum_undir_dist"]]) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_origingridnodatawindows")
}
