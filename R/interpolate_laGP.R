#' interpolate_laGP
#'
#' @param independent test
#' @param dependent test
#' @param pred_grid test
#' @param auto test
#' @param d test
#' @param g test
#' @param on_residuals test
#'
#' @return test
#'
#' @examples
#' \donttest{
#' independent <- tibble::tribble(
#'   ~x, ~y, ~z,
#'   1,1,2,
#'   3,2,1,
#'   2,3,3,
#'   7,7,7,
#'   8,7,9,
#'   9,7,8
#' )
#'
#' dependent <- c(1,2,1,5,5,6)
#'
#' pred_grid <- tibble::as_tibble(expand.grid(x = 1:10, y = 1:10, z = 1:10))
#'
#' pred <- interpolate_laGP(independent, dependent, pred_grid, auto = F, d = c(3, 3, 4), g = 0.1, on_residuals = F)
#'
#' pred_grid$pred_mean <- pred$mean
#'
#' ggplot(data = pred_grid[pred_grid$z == 4, ]) +
#'   geom_raster(aes(x, y, fill = pred_mean))
#'
#' }
#'
#' @export
interpolate_laGP <- function(independent, dependent, pred_grid, auto = T, d, g, on_residuals = F) {

  if (on_residuals) {
    # linear fit
    model <- stats::lm(dependent ~ x + y + z, data = independent)
    dependent <- model[["residuals"]]
  }

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

  if (on_residuals) {
    # linear fit
    pred$mean <- pred$mean + stats::predict(model, pred_grid[c("x", "y", "z")])
  }

  # return result
  return(pred)
}

# rescaling helpers
range_01 <- function(x, min, max) { (x - min) / (max - min) }
dist_scale_01 <- function(x, min, max) { x / abs(min - max) }
