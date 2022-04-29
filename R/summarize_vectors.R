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
#' @rdname average_origin
NULL

#' @rdname average_origin
#' @export
average_origin_searchid <- function(origin_grid) {
  UseMethod("average_origin_searchid")
}

#' @rdname average_origin
#' @export
average_origin_searchid.default <- function(origin_grid) {
  stop("x is not an object of class mobest_origingrid")
}

#' @rdname average_origin
#' @export
average_origin_searchid.mobest_origingrid <- function(origin_grid) {
  origin_grid %>%
    dplyr::group_by(.data[["search_id"]]) %>%
    dplyr::summarise(
      mean_search_z = mean(.data[["search_z"]]),
      sd_search_z = stats::sd(.data[["search_z"]]),
      region_id = dplyr::first(.data[["region_id"]]),
      undirected_mean_spatial_distance = mean(.data[["spatial_distance"]]),
      undirected_sd_spatial_distance = stats::sd(.data[["spatial_distance"]]),
      directed_mean_spatial_distance = sqrt(
        mean(.data[["search_x"]] - .data[["origin_x"]])^2 +
          mean(.data[["search_y"]] - .data[["origin_y"]])^2
      ) / 1000,
      mean_angle_deg = mobest::vec2deg(
        c(mean(.data[["origin_x"]] - .data[["search_x"]]), mean(.data[["origin_y"]] - .data[["search_y"]]))
      ),
      mean_angle_deg_cut = cut_angle_deg(.data[["mean_angle_deg"]]),
      .groups = "drop"
    ) %>% dplyr::arrange(.data[["undirected_mean_spatial_distance"]]) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_meanorigingrid")
}

#' @rdname average_origin
#' @export
average_origin_moving_window <- function(
  origin_grid, window_start, window_stop, window_width, window_step
) {
  UseMethod("average_origin_moving_window")
}

#' @rdname average_origin
#' @export
average_origin_moving_window.default <- function(
  origin_grid, window_start, window_stop, window_width, window_step
) {
  stop("x is not an object of class mobest_origingrid")
}

#' @rdname average_origin
#' @export
average_origin_moving_window.mobest_origingrid <- function(
  origin_grid, window_start, window_stop, window_width, window_step
) {
  future::plan(future::multisession)
  furrr::future_map_dfr(
    # loop through all regions
    unique(origin_grid$region_id),
    function(region) {
      origin_per_region <- origin_grid %>%
        dplyr::filter(.data[["region_id"]] == region)
      purrr::map2_df(
        # define moving windows and loop through them
        seq(window_start, window_stop - window_width, window_step),
        seq(window_start + window_width, window_stop, window_step),
        function(start, end) {
          # prepare window data subsets
          io <- dplyr::filter(
            origin_per_region,
            .data[["search_z"]] >= start,
            .data[["search_z"]] < end
          )
          io_upper_quartile <- dplyr::filter(
            io,
            .data[["spatial_distance"]] >= stats::quantile(.data[["spatial_distance"]], probs = 0.75)
          )
          io_run_grouped <- io %>%
            dplyr::group_by(.data[["search_id"]]) %>%
            dplyr::summarise(
              mean_spatial_distance = mean(.data[["spatial_distance"]]),
              groups = "drop"
            )
          if (nrow(io) > 0) {
            tibble::tibble(
              z = mean(c(start, end)),
              region_id = region,
              undirected_mean_spatial_distance =
                mean(io$spatial_distance),
              undirected_mean_spatial_distance_upper_quartile =
                mean(io_upper_quartile$spatial_distance),
              directed_mean_spatial_distance = sqrt(
                mean(io$search_x - io$origin_x)^2 +
                  mean(io$search_y - io$origin_y)^2
              ) / 1000,
              directed_mean_spatial_distance_upper_quartile = sqrt(
                mean(io_upper_quartile$search_x - io_upper_quartile$origin_x)^2 +
                  mean(io_upper_quartile$search_y - io_upper_quartile$origin_y)^2
              ) / 1000,
              se_spatial_distance = if (nrow(io_run_grouped) >= 3) {
                se(io_run_grouped$mean_spatial_distance)
              } else {
                Inf
              },
              sd_spatial_distance = if (nrow(io_run_grouped) >= 3) {
                stats::sd(io_run_grouped$mean_spatial_distance)
              } else {
                Inf
              },
              fraction_smaller_500 = sum(io$spatial_distance < 500) / nrow(io),
              fraction_bigger_500 = sum(io$spatial_distance >= 500 & io$spatial_distance < 1000) / nrow(io),
              fraction_bigger_1000 = sum(io$spatial_distance >= 1000 & io$spatial_distance < 2000) / nrow(io),
              fraction_bigger_2000 = sum(io$spatial_distance >= 2000) / nrow(io)
            )
          } else {
            tibble::tibble(
              z = mean(c(start, end)),
              region_id = region,
              undirected_mean_spatial_distance = NA,
              directed_mean_spatial_distance = NA,
              mean_angle_deg = NA,
              se_spatial_distance = Inf,
              sd_spatial_distance = Inf,
              fraction_smaller_500 = NA,
              fraction_bigger_500 = NA,
              fraction_bigger_1000 = NA,
              fraction_bigger_2000 = NA
            )
          }
        }
      )
    }
  ) %>%
  tibble::new_tibble(., nrow = nrow(.), class = "mobest_movingorigingrid")
}

se <- function(x) stats::sd(x)/sqrt(length(x))

#' @rdname average_origin
#' @export
no_data_windows <- function(moving_origin_grid, window_step) {
  UseMethod("no_data_windows")
}

#' @rdname average_origin
#' @export
no_data_windows.default <- function(moving_origin_grid, window_step) {
  stop("x is not an object of class mobest_movingorigingrid")
}

#' @rdname average_origin
#' @export
no_data_windows.mobest_movingorigingrid <- function(moving_origin_grid, window_step) {
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
