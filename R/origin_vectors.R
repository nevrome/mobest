#' Search for a past point of "origin" in a prediction grid
#'
#' @param locate_product An object of class \code{mobest_locateproduct} as created by
#' \link{locate} + \link{locate_multi} and \link{multiply_dependent_probabilities}.
#' @param ... (Additional) grouping variables (\code{independent_table_id}, \code{dependent_setting_id},
#' \code{kernel_setting_id}, \code{pred_grid_id}, ...)
#' @param quiet Logical. Should a progress indication be printed?
#'
#' @return An object of class \code{mobest_origin_grid}
#'
#' @rdname origin_vectors
#' @export
determine_origin_vectors <- function(
  locate_product,
  ...,
  quiet = F
) {
  .grouping_var <- rlang::ensyms(...)
  # input checks
  checkmate::assert_class(locate_product, "mobest_locateproduct")
  # summarise data
  locate_groups <- locate_product %>%
    dplyr::group_split(
      !!!.grouping_var,
      .data[["search_id"]]
    )
  origin_grid <- locate_groups %>%
    purrr::map_df(
      function(locate_group) {
      locate_group %>%
        dplyr::mutate(
          ov_x = .data[["field_x"]] - .data[["search_x"]],
          ov_y = .data[["field_y"]] - .data[["search_y"]],
          ov_dist = sqrt(.data[["ov_x"]]^2 + .data[["ov_y"]]^2),
          ov_dist_sd = sqrt(wtd.var(.data[["ov_dist"]], .data[["probability"]]))
        ) %>%
        # keep only the maximum probability grid cell
        dplyr::slice_max(.data[["probability"]], n = 1, with_ties = FALSE) %>%
        dplyr::mutate(
          ov_angle_deg = vec2deg(c(.data[["ov_x"]], .data[["ov_y"]])),
          ov_angle_cut = cut_angle_deg(.data[["ov_angle_deg"]])
        )
      }
    )
  # compile output
  origin_grid %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_originvectors")
}

#### helper functions ####

# This function was copied from the Hmisc R package, v4.7-0
# https://CRAN.R-project.org/package=Hmisc
wtd.var <- function (
  x, weights = NULL, normwt = FALSE, na.rm = TRUE,
  method = c("unbiased", "ML")
) {
  method <- match.arg(method)
  if (!length(weights)) {
    if (na.rm)
      x <- x[!is.na(x)]
    return(stats::var(x))
  }
  if (na.rm) {
    s <- !is.na(x + weights)
    x <- x[s]
    weights <- weights[s]
  }
  if (normwt)
    weights <- weights * length(x)/sum(weights)
  if (normwt || method == "ML")
    return(as.numeric(stats::cov.wt(cbind(x), weights, method = method)$cov))
  sw <- sum(weights)
  if (sw <= 1)
    warning("only one effective observation; variance estimate undefined")
  xbar <- sum(weights * x)/sw
  sum(weights * ((x - xbar)^2))/(sw - 1)
}
