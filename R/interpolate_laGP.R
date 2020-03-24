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
  # priors for the global GP
  if (auto) {
    da <- laGP::darg(list(mle = TRUE, max=10), independent)
    ga <- laGP::garg(list(mle = TRUE, max=10), dependent)
    d <- da$start
    g <- ga$start
  }
  # fit the global GP
  gp <- laGP::newGPsep(X = independent, Z = dependent, d = d, g = g, dK = auto)
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
  pred <- laGP::predGPsep(gp, XX = pred_grid[, c("x_01", "y_01", "z_01")], lite = T)
  # delete GP object
  laGP::deleteGPsep(gp)
  # return result
  return(pred)
}
