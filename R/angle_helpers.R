#' Angle transformation functions
#'
#' Functions to transform angles from and to different formats.
#' The functions are not vectorised.
#'
#' @param x Input for transformation
#' \itemize{
#'  \item{deg: }{Double. Angle in degrees}
#'  \item{rad: }{Double. Angle in radians}
#'  \item{vec: }{2 element double vector. Angle as unit vector}
#' }
#'
#' @return Same formats as x
#'
#' @rdname angle_transformers
#' @export
deg2rad <- function(x) {
  x * pi/180
}

#' @rdname angle_transformers
#' @export
rad2deg <- function(x) {
  x * 180/pi
}

#' @rdname angle_transformers
#' @export
deg2vec <- function(x) {
  c(sin(deg2rad(x)), cos(deg2rad(x)))
}

#' @rdname angle_transformers
#' @export
vec2deg <- function(x) {
  res <- rad2deg(atan2(x[1], x[2]))
  if (res < 0) {
    360 + res
  } else {
    res
  }
}

#' @rdname angle_transformers
#' @export
rad2vec <- function(x) {
  c(sin(x), cos(x))
}

#' @rdname angle_transformers
#' @export
vec2rad <- function(x) {
  atan2(x[1], x[2])
}

#' Mean angle functions
#'
#' Functions to calculate the mean of multiple angles.
#'
#' @param x Double vector. Angles in degrees
#'
#' @return 2 element double vector. Angle as unit vector
#'
#' @rdname angle_helpers
#' @export
mean_deg2vec <- function(x) {
  y <- lapply(x, deg2vec)
  Reduce(`+`, y)/length(y)
}
