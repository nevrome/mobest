#' interpolate_laGP
#'
#' @param independent test
#' @param dependent test
#' @param pred_grid test
#' @param auto test
#' @param d test
#' @param g test
#'
#' @return test
#'
#' @export
interpolate_laGP <- function(independent, dependent, pred_grid, auto = T, d, g) {

  minx <- min(independent[["x"]])
  maxx <- max(independent[["x"]])
  miny <- min(independent[["y"]])
  maxy <- max(independent[["y"]])
  minz <- min(independent[["z"]])
  maxz <- max(independent[["z"]])

  # rescale
  independent_rescaled <- independent %>%
    dplyr::mutate(
      x = range_01(.data[["x"]], minx, maxx),
      y = range_01(.data[["y"]], miny, maxy),
      z = range_01(.data[["z"]], minz, maxz)
    )

  pred_grid_rescaled <- pred_grid %>%
    dplyr::mutate(
      x = range_01(.data[["x"]], minx, maxx),
      y = range_01(.data[["y"]], miny, maxy),
      z = range_01(.data[["z"]], minz, maxz)
    )

  d_rescaled <- d
  d_rescaled[1] <- dist_scale_01(d[1], minx, maxx)
  d_rescaled[2] <- dist_scale_01(d[2], miny, maxy)
  d_rescaled[3] <- dist_scale_01(d[3], minz, maxz)

  # priors for the global GP
  if (auto) {
    da <- laGP::darg(list(mle = TRUE, max = 10), independent_rescaled)
    ga <- laGP::garg(list(mle = TRUE, max = 10), dependent)
    d_rescaled <- da$start
    g <- ga$start
  }

  # fit the global GP
  gp <- laGP::newGPsep(X = independent_rescaled, Z = dependent, d = d_rescaled, g = g, dK = auto)

  # optimise fit automatically
  if (auto) {
    laGP::mleGPsep(
      gpsepi = gp,
      param = "both",
      tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab),
      maxit = 200
    )
  }

  # predictions from the global GP on the prediction
  pred <- laGP::predGPsep(gp, XX = pred_grid_rescaled[, c("x", "y", "z")], lite = T)

  # delete GP object
  laGP::deleteGPsep(gp)

  # return result
  return(pred)
}

# rescaling helpers
range_01 <- function(x, min, max) { (x - min) / (max - min) }
dist_scale_01 <- function(x, min, max) { x / abs(min - max) }
