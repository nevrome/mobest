#' Use laGPs mle algorithms for kriging kernel parameter estimation
#'
#' @param independent An object of class mobest_spatiotemporalpositions
#' @param dependent An object of class mobest_observations
#' @param iterations Integer. Number of mle iterations
#' @param total_scaling_factor Numeric. It turned out that the parameter estimation
#' works better in certain numeric ranges. This scaling factor affects all three
#' dimensions
#' @param g Numeric. Fixed value for the nugget term
#' @param space_time_scaling_factor_sequence Numeric vector. Sequence of space-time
#' scaling factors
#' @param verb Integer. See \link[laGP]{mleGP} for more information
#'
#' @rdname parameter_estimation_mle
#' @export
laGP_mle_sequence_isotropic_fixed_g <- function(
  independent, dependent, iterations, total_scaling_factor = 1000,
  g, space_time_scaling_factor_sequence, verb = 1) {
  # input check and transformation
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  independent_without_id <- independent %>%
    dplyr::transmute(
      x = x / total_scaling_factor,
      y = y / total_scaling_factor,
      z = z / total_scaling_factor
    )
  checkmate::assert_class(dependent, "mobest_observations")
  checkmate::assert_count(iterations)
  checkmate::assert_count(total_scaling_factor)
  checkmate::assert_numeric(space_time_scaling_factor_sequence)
  # run parameter estimation loop for each dependent variable and iteration
  purrr::map2_dfr(
    dependent, names(dependent),
    function(cur_dependent, cur_dependent_name) {
    # run multiple iterations
    purrr::map_dfr(1:iterations, function(i) {
      # run for each scaling factor
      mleGP_out_list <- purrr::map(space_time_scaling_factor_sequence, function(scaling_factor) {
        independent_rescaled <- independent_without_id %>%
          dplyr::mutate(
            z = z * scaling_factor
          )
        # parameter estimation
        da <- laGP::darg(list(mle = TRUE), independent_rescaled)
        gp <- laGP::newGP(
          X = independent_rescaled,
          Z = cur_dependent,
          d = da$start,
          g = g,
          dK = TRUE
        )
        param_estimation <- laGP::mleGP(
          gpi = gp,
          param = "d",
          tmin = da$min, tmax = da$max, ab = da$ab,
          verb = verb
        )
        param_estimation$l <- laGP::llikGP(gp, dab = da$ab)
        laGP::deleteGP(gp)
        return(param_estimation)
      })
      # combine list to better readable data.frame
      tibble::tibble(
        iteration = i,
        ancestry_component = cur_dependent_name,
        scaling_factor = space_time_scaling_factor_sequence,
        scaling_factor_fractional = fractional::fractional(space_time_scaling_factor_sequence),
        scaling_factor_label = factor(
          as.character(as.character(scaling_factor_fractional)),
          levels = as.character(as.character(scaling_factor_fractional))
        ),
        d = sqrt(sapply(mleGP_out_list, function(x) { x$d })) * total_scaling_factor,
        l = sapply(mleGP_out_list, function(x) { x$l }),
        its = sapply(mleGP_out_list, function(x) { x$it })
      ) %>% dplyr::mutate(
        ds = d,
        dt = d / scaling_factor
      )
    })
  })
}

#' @rdname parameter_estimation_mle
#' @export
laGP_mle_anisotropic <- function(
  independent, dependent, iterations, total_scaling_factor = 1000, verb = 1) {
  # input check and transformation
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  independent_without_id <- independent %>%
    dplyr::transmute(
      x = x / total_scaling_factor,
      y = y / total_scaling_factor,
      z = z / total_scaling_factor
    )
  checkmate::assert_class(dependent, "mobest_observations")
  checkmate::assert_count(iterations)
  checkmate::assert_count(total_scaling_factor)
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
        maxit = 200,
        verb = verb
      )
      laGP::deleteGPsep(gp)
      return(param_estimation)
    })
    # look at result parameters
    purrr::map_dfr(mleGPsep_params, function(x) {
      tibble::tibble(
        mle_method = "mleGPsep",
        ancestry_component = cur_dependent_name,
        dx = sqrt(x$theta[1]) * total_scaling_factor,
        dy = sqrt(x$theta[2]) * total_scaling_factor,
        dt = sqrt(x$theta[3]) * total_scaling_factor,
        g = x$theta[4],
        its = x$its,
        msg = x$msg,
        conv = x$conv
      )
    })
  })
}

#' @rdname parameter_estimation_mle
#' @export
laGP_jmle_anisotropic <- function(
  independent, dependent, iterations, total_scaling_factor = 1000, verb = 1) {
  # input check and transformation
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  independent_without_id <- independent %>%
    dplyr::transmute(
      x = x / total_scaling_factor,
      y = y / total_scaling_factor,
      z = z / total_scaling_factor
    )
  checkmate::assert_class(dependent, "mobest_observations")
  checkmate::assert_count(iterations)
  checkmate::assert_count(total_scaling_factor)
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
          maxit = 200,
          verb = verb
        )
        laGP::deleteGPsep(gp)
        return(param_estimation)
      })
      # look at result parameters
      jmleGPsep_params %>%
        dplyr::transmute(
          mle_method = "jmleGPsep",
          ancestry_component = cur_dependent_name,
          dx = sqrt(d.1) * total_scaling_factor,
          dy = sqrt(d.2) * total_scaling_factor,
          dt = sqrt(d.3) * total_scaling_factor,
          g = g,
          its = tot.its,
          msg = NA,
          conv = dconv
        ) %>%
        tibble::as_tibble()
    })
}
