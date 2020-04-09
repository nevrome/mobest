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
#'   0,0,0,
#'   1,0,0,
#'   0,1,0,
#'   0,0,1,
#'   10,10,10,
#'   9,10,10,
#'   10,9,10,
#'   10,10,9
#' )
#'
#' dependent <- c(1,1,1,1,10,10,10,10)
#'
#' pred_grid <- tibble::as_tibble(expand.grid(x = 0:10, y = 0:10, z = 0:10))
#'
#' pred <- interpolate_laGP(independent, dependent, pred_grid, auto = F, d = c(3, 3, 3)^2, g = 0.01, on_residuals = T)
#'
#' pred_grid$pred_mean <- pred$mean
#'
#' library(ggplot2)
#' ggplot(data = pred_grid) +
#'   geom_raster(aes(x, y, fill = pred_mean)) +
#'   facet_wrap(~z) +
#'   scale_fill_viridis_c()
#'
#' }
#'
#' @export
interpolate_laGP <- function(independent, dependent, pred_grid, auto = F, d, g, on_residuals = T) {

  if (on_residuals) {
    # linear fit
    combined <- independent %>% dplyr::mutate(d = dependent)
    model <- stats::lm(d ~ x + y + z, data = combined)
    dependent <- model[["residuals"]]
  }

  # priors for the global GP
  if (auto) {
    da <- laGP::darg(list(mle = TRUE, max = 10), independent)
    ga <- laGP::garg(list(mle = TRUE, max = 10), dependent)
    d <- da$start
    g <- ga$start
  }

  # fit the global GP
  gp <- laGP::newGPsep(X = independent, Z = dependent, d = d, g = g)

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
  pred <- laGP::predGPsep(gp, XX = pred_grid[, c("x", "y", "z")], lite = T)

  # delete GP object
  laGP::deleteGPsep(gp)

  if (on_residuals) {
    # add predictions from linear model again
    pred$mean <- pred$mean + stats::predict(model, pred_grid[c("x", "y", "z")])
  }

  # return result
  return(pred)
}
