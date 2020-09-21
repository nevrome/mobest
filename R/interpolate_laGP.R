#' 3D interpolation with laGP
#'
#' Performs kriging with laGPs gaussian process prediction. See
#' \code{?laGP::predGPsep} for more information. If \code{on_residuals = T} a linear model
#' is wrapped around the kriging model to handle the main trends independently.
#'
#' @param independent Dataframe with input point position coordinates x, y and z
#' @param dependent Vector with input point values
#' @param pred_grid Dataframe with output point position coordinates x, y and z
#' @param d Numeric vector. Lengthscale parameter. See \code{?laGP::newGP} for more info
#' @param g Numeric. Nugget parameter
#' @param auto Should the lengthscale and nugget values be automatically determined
#' by laGPs maximum likelihood algorithm? See \code{?laGP::mleGPsep} for more info
#' @param on_residuals Should a linear model take out the main trends before the kriging interpolation?
#'
#' @return Output of \code{?laGP::predGPsep}
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
interpolate_laGP <- function(independent, dependent, pred_grid, d, g, auto = F, on_residuals = T) {

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
