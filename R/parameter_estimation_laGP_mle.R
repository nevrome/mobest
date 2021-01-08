#' Title
#'
#' @param independent
#' @param dependent
#' @param iterations
#' @param scaling_factor
#'
#' @return
#' @export
laGP_mle_anisotropic <- function(independent, dependent, iterations, scaling_factor = 1000) {
  # input check and transformation
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  independent_without_id <- independent %>%
    dplyr::transmute(
      x = x / scaling_factor,
      y = y / scaling_factor,
      z = z / scaling_factor
    )
  checkmate::assert_class(dependent, "mobest_observations")
  checkmate::assert_count(iterations)
  checkmate::assert_count(scaling_factor)
  # run parameter estimation loop for each dependent variable and iteration
  purrr::map2_dfr(
    dependent, names(dependent),
    function(cur_dependent, cur_dependent_name) {
    # parameter estimation
    mleGPsep_params <- purrr::map(1:iterations, function(i) {
      da <- laGP::darg(list(mle = TRUE), independent_without_id)
      ga <- laGP::garg(list(mle = TRUE), cur_dependent)
      gp <- laGP::newGPsep(
        X = independent_without_id,
        Z = cur_dependent,
        d = da$start,
        g = ga$start,
        dK = TRUE
      )
      param_estimation <- laGP::mleGPsep(
        gpsepi = gp,
        param = "both",
        tmin = c(da$min, ga$min), tmax = c(da$max, ga$max), ab = c(da$ab, ga$ab),
        maxit = 200
      )
      laGP::deleteGPsep(gp)
      return(param_estimation)
    })
    # look at result parameters
    purrr::map_dfr(mleGPsep_params, function(x) {
      tibble::tibble(
        mle_method = "mleGPsep",
        ancestry_component = cur_dependent_name,
        dx = sqrt(x$theta[1]) * scaling_factor,
        dy = sqrt(x$theta[2]) * scaling_factor,
        dt = sqrt(x$theta[3]) * scaling_factor,
        g = x$theta[4],
        its = x$its,
        msg = x$msg,
        conv = x$conv
      )
    })
  })
}

#' Title
#'
#' @param independent
#' @param dependent
#' @param iterations
#' @param scaling_factor
#'
#' @return
#' @export
laGP_jmle_anisotropic <- function(independent, dependent, iterations, scaling_factor = 1000) {
  # input check and transformation
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  independent_without_id <- independent %>%
    dplyr::transmute(
      x = x / scaling_factor,
      y = y / scaling_factor,
      z = z / scaling_factor
    )
  checkmate::assert_class(dependent, "mobest_observations")
  checkmate::assert_count(iterations)
  checkmate::assert_count(scaling_factor)
  # run parameter estimation loop for each dependent variable and iteration
  purrr::map2_dfr(
    dependent, names(dependent),
    function(cur_dependent, cur_dependent_name) {
      # parameter estimation
      jmleGPsep_params <- purrr::map_dfr(1:iterations, function(i) {
        da <- laGP::darg(list(mle = TRUE), independent_without_id)
        ga <- laGP::garg(list(mle = TRUE), cur_dependent)
        gp <- laGP::newGPsep(
          X = independent_without_id,
          Z = cur_dependent,
          d = da$start,
          g = ga$start,
          dK = TRUE
        )
        param_estimation <- laGP::jmleGPsep(
          gpsepi = gp,
          drange = c(da$min, da$max),
          grange = c(ga$min, ga$max),
          dab = da$ab,
          gab = ga$ab,
          maxit = 200
        )
        laGP::deleteGPsep(gp)
        return(param_estimation)
      })
      # look at result parameters
      jmleGPsep_params %>%
        dplyr::transmute(
          mle_method = "jmleGPsep",
          ancestry_component = cur_dependent_name,
          dx = sqrt(d.1) * scaling_factor,
          dy = sqrt(d.2) * scaling_factor,
          dt = sqrt(d.3) * scaling_factor,
          g = g,
          its = tot.its,
          msg = NA,
          conv = dconv
        ) %>%
        tibble::as_tibble()
    })
}
