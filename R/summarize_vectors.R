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
#' @param window_start Numeric. Start date of the moving window sequence
#' @param window_stop Numeric. Stop date of the moving window sequence
#' @param window_width Numeric. Width of each individual moving window
#' @param window_step Numeric. Frequency of moving windows. Example:
#' If the first window starts at -3500 and extends until -3200, should the next
#' one start at -3450 or -3400, so with window_step = 50 or window_step = 100?
#' @param dist_fraction_breaks Numeric vector. Cutting breaks of distance fractions to cross-tabulate
#' for each moving window. If set, then a list column \code{ov_dist_fractions} is added to the output,
#' which features a tibble with the sample counts and fractions for each moving window
#' \code{ov_dist_fractions}
#' @param quiet Logical. Should a progress indication be printed?
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
  origin_vectors, ..., window_start, window_stop, window_width, window_step,
  dist_fraction_breaks = NA,
  quiet = F
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
    dplyr::group_split(!!!.grouping_var)
  # loop through units
  if (!quiet) { message("Summarising groups") }
  origin_summary <- purrr::imap_dfr(
    vector_groups,
    function(vector_group, i) {
      # define moving windows and loop through them
      windows_starts <- seq(window_start, window_stop - window_width, window_step)
      windows_stops  <- seq(window_start + window_width, window_stop, window_step)
      if (!quiet) {
        message("Summarising windows for group ", i, " of ", length(vector_groups))
        pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = length(windows_starts))
      }
      purrr::map2_df(
        windows_starts, windows_stops,
        function(start, end) {
          if (!quiet) { pb$tick() }
          # prepare window data subsets
          io <- dplyr::filter(
            vector_group,
            .data[["search_z"]] >= start,
            .data[["search_z"]] < end
          )
          # check if there are any samples in this window
          if (nrow(io) > 0) {
            io %>%
              dplyr::group_by(!!!.grouping_var) %>%
              dplyr::summarise(
                z            = mean(c(start, end)),
                ov_dist      = sqrt(mean(.data[["ov_x"]])^2 + mean(.data[["ov_y"]])^2),
                ov_dist_se   = if (dplyr::n() >= 2) {
                                 calculate_standard_error(sqrt(.data[["ov_x"]]^2 + .data[["ov_y"]]^2))
                               } else { Inf },
                ov_dist_sd   = if (dplyr::n() >= 2) {
                                 stats::sd(sqrt(.data[["ov_x"]]^2 + .data[["ov_y"]]^2))
                               } else { Inf },
                ov_angle_deg = vec2deg(c(mean(.data[["ov_x"]]), mean(.data[["ov_x"]]))),
                ov_dist_fractions =
                  # check if dist_fraction_breaks is set
                  if (!any(is.na(dist_fraction_breaks))) {
                    list(determine_dist_fractions(sqrt(.data[["ov_x"]]^2 + .data[["ov_y"]]^2), dist_fraction_breaks))
                  } else {
                    NULL
                  }
              )
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
                ov_angle_deg = NA,
                ov_dist_fractions =
                  if (!any(is.na(dist_fraction_breaks))) {
                    # to get a proper df back for every time window
                    list(determine_dist_fractions(c(), dist_fraction_breaks))
                  } else {
                    NULL
                  }
              )
          }
        }
      )
    }
  )
  if (!quiet) { message("Compiling output") }
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

determine_dist_fractions <- function(x, cutting_points) {
  x_cut <- if (length(x) != 0) {
    cutting_points[cut(x, cutting_points, labels = F, include.lowest = T)]
  } else {
    c()
  }
  # to get all fraction levels for each window, even if they don't occur
  x_cut_factor <- factor(x_cut, levels = head(cutting_points, -1))
  x_cut_factor %>%
    base::table() %>%
    unclass() %>%
    tibble::enframe(name = "ov_dist_length_lower_end", value = "count") %>%
    dplyr::mutate(
      ov_dist_length_lower_end = as.numeric(.data[["ov_dist_length_lower_end"]]),
      fraction = if (sum(.data[["count"]]) != 0) {
        .data[["count"]]/sum(.data[["count"]])
      } else {
        0
      }
    ) %>%
    tibble::add_column(
      ov_dist_length_upper_end = tail(cutting_points, -1),
      .after = "ov_dist_length_lower_end"
    )
}
