#' angle helpers
#' @name angle_helpers
#'
#' @param x test
#'
#' @rdname angle_helpers
#'
#' @export
deg2rad <- function(x) {
  x * pi/180
}

#' @rdname angle_helpers
#' @export
rad2deg <- function(x) {
  x * 180/pi
}

#' @rdname angle_helpers
#' @export
deg2vec <- function(x) {
  c(sin(deg2rad(x)), cos(deg2rad(x)))
}

#' @rdname angle_helpers
#' @export
vec2deg <- function(x) {
  res <- rad2deg(atan2(x[1], x[2]))
  if (res < 0) {
    360 + res
  } else {
    res
  }
}

#' @rdname angle_helpers
#' @export
mean_vec <- function(x) {
  y <- lapply(x, deg2vec)
  Reduce(`+`, y)/length(y)
}

#' @rdname angle_helpers
#' @export
mean_deg <- function(x) {
  vec2deg(mean_vec(x))
}
