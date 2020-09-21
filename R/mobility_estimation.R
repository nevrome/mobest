#' estimate_mobility
#'
#' @param origin_grid Either an object of class \code{mobest_origin_grid} as
#' produced by \link{search_spatial_origin} or just an object of class
#' \code{mobest_interpol_grid}. Depending on the input a different mobility
#' estimation algorithm is applied
#' @param mobility_regions test
#'
#' @return test
#'
#' @rdname estimate_mobility
#' @export
estimate_mobility <- function(origin_grid, mobility_regions) {
  UseMethod("estimate_mobility")
}

#' @rdname estimate_mobility
#' @export
estimate_mobility.default <- function(origin_grid, mobility_regions) {
  stop("x is not an object of class mobest_interpol_grid or mobest_origin_grid")
}

#' @rdname estimate_mobility
#' @export
estimate_mobility.mobest_origin_grid <- function(origin_grid, mobility_regions) {



}
