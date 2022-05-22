#' Create and run a kriging model grid
#'
#' Construct and run a model grid with all combinations of the different input parameter
#' configurations for the kriging model.
#'
#' @param independent An object of class \code{mobest_spatiotemporalpositions_multi}
#' as created by \link{create_spatpos_multi}. Spatiotemporal input point positions.
#' @param dependent An object of class \code{mobest_observations_multi}
#' as created by \link{create_obs_multi}. Dependent variables that should be interpolated.
#' @param kernel An object of class \code{mobest_kernelsetting_multi}
#' as created by \link{create_kernset_multi}. Kernel parameter settings.
#' @param prediction_grid An object of class \code{mobest_spatiotemporalpositions_multi}
#' as for example created by \link{create_prediction_grid}.
#' Prediction grid positions for the interpolation.
#' @param model_grid An object of class \code{mobest_modelgrid}
#' as created by \link{create_model_grid}.
#' @param unnest Logical. Should the kriging result be unnested to return a
#' prediction point-wise table of class \code{mobest_interpolgrid}?
#' @param quiet Logical. Should a progress indication be printed?
#'
#' @name interpolate
NULL

#' @rdname interpolate
#' @export
create_model_grid <- function(
  independent,
  dependent,
  kernel,
  prediction_grid
) {
  # input check
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions_multi")
  checkmate::assert_class(dependent, "mobest_observations_multi")
  checkmate::assert_class(kernel, "mobest_kernelsetting_multi")
  checkmate::assert_class(prediction_grid, "mobest_spatiotemporalpositions_multi")
  check_compatible_multi(kernel, dependent, check_names_equal)
  check_compatible_multi(independent, dependent, check_df_nrow_equal)
  # fill create general structure and id columns
  independent_tables <- tibble::tibble(
    independent_table = independent,
    independent_table_id = names(independent)
  )
  dependent_vars <- purrr::map2_dfr(
    names(dependent), dependent,
    function(dependent_name, one_dependent) {
      tibble::tibble(
        dependent_setting_id = dependent_name,
        dependent_var_id = names(one_dependent),
        dependent_var = as.list(one_dependent)
      )
    }
  )
  kernel_settings <- purrr::map2_dfr(
    names(kernel), kernel,
    function(kernel_name, one_kernel) {
      tibble::tibble(
        kernel_setting_id = kernel_name,
        dependent_var_id = names(one_kernel),
        kernel_setting = one_kernel[1:length(one_kernel)]
      )
    }
  )
  pred_grids <- tibble::tibble(
    pred_grid = prediction_grid,
    pred_grid_id = names(prediction_grid)
  )
  # expand grid and create model grid
  model_grid <- create_model_grid_raw(
    independent_tables = independent_tables,
    dependent_vars = dependent_vars,
    kernel_settings = kernel_settings,
    pred_grids = pred_grids
  ) %>%
    # make subclass of tibble
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_modelgrid")
  return(model_grid)
}

create_model_grid_raw <- function(independent_tables, dependent_vars, kernel_settings, pred_grids) {
  expand.grid(
    independent_table_id = unique(independent_tables[["independent_table_id"]]),
    dependent_setting_id = unique(dependent_vars[["dependent_setting_id"]]),
    dependent_var_id = unique(dependent_vars[["dependent_var_id"]]),
    kernel_setting_id = unique(kernel_settings[["kernel_setting_id"]]),
    pred_grid_id = unique(pred_grids[["pred_grid_id"]]),
    stringsAsFactors = F
  ) %>%
    dplyr::left_join(
      independent_tables, by = "independent_table_id"
    ) %>%
    dplyr::left_join(
      dependent_vars, by = c("dependent_setting_id", "dependent_var_id")
    ) %>%
    dplyr::left_join(
      kernel_settings, by = c("kernel_setting_id", "dependent_var_id")
    ) %>%
    dplyr::left_join(
      pred_grids, by = "pred_grid_id"
    )
}

#' @rdname interpolate
#' @export
run_model_grid <- function(model_grid, unnest = T, quiet = F) {
  # input check
  checkmate::assert_class(model_grid, "mobest_modelgrid")
  # run interpolation for each entry in the model_grid
  if (!quiet) {
    message("Running models")
    pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = nrow(model_grid))
  }
  prediction <- purrr::map(
    1:nrow(model_grid), function(i) {
      if (!quiet) { pb$tick() }
      interpolate(
        independent = model_grid[["independent_table"]][[i]],
        dependent = model_grid[["dependent_var"]][[i]],
        pred_grid = model_grid[["pred_grid"]][[i]],
        # d has to be squared because of the configuration of the default laGP kernel
        d = as.numeric(model_grid[["kernel_setting"]][[i]][c("dsx", "dsy", "dt")])^2,
        g = model_grid[["kernel_setting"]][[i]][["g"]],
        auto = model_grid[["kernel_setting"]][[i]][["auto"]],
        on_residuals = model_grid[["kernel_setting"]][[i]][["on_residuals"]]
      )
    }
  )
  # simplify model_grid
  if (!quiet) { message("Finishing prediction output") }
  model_grid_simplified <- model_grid %>%
    dplyr::mutate(
      dsx = purrr::map_dbl(model_grid$kernel_setting, purrr::pluck("dsx")),
      dsy = purrr::map_dbl(model_grid$kernel_setting, purrr::pluck("dsy")),
      dt = purrr::map_dbl(model_grid$kernel_setting, purrr::pluck("dt")),
      g = purrr::map_dbl(model_grid$kernel_setting, purrr::pluck("g"))
    ) %>%
    dplyr::select(
      -.data[["kernel_setting"]],
      -.data[["independent_table"]],
      -.data[["dependent_var"]],
      -.data[["pred_grid"]]
    )
  # add prediction results for each run as a data.frame in a list column to model_grid
  model_grid_simplified$prediction <- purrr::map2(
    prediction, model_grid[["pred_grid"]],
    function(x, y) {
      pred <- data.frame(
        mean = x$mean,
        sd = sqrt(x$s2),
        stringsAsFactors = F
      )
      # merge with pred_grid to add relevant spatial information
      cbind(y, pred)
    }
  )
  if (unnest) {
    model_grid_simplified %>%
      tidyr::unnest(cols = "prediction") %>%
      # make subclass of tibble
      tibble::new_tibble(., nrow = nrow(.), class = "mobest_interpolgrid") %>%
      return()
  } else {
    return(model_grid_simplified)
  }
}

#' 3D interpolation with laGP
#'
#' Performs kriging with laGPs gaussian process prediction. See
#' \code{?laGP::predGPsep} for more information.
#'
#' @param independent An object of class mobest_spatiotemporalpositions.
#' Spatiotemporal input point positions
#' @param dependent Numeric vector.
#' Dependent variable that should be interpolated
#' @param pred_grid An object of class mobest_spatiotemporalpositions.
#' Prediction grid positions for the interpolation
#' @param d Numeric vector. Lengthscale parameter.
#' See \code{?laGP::newGP} for more info
#' @param g Numeric. Nugget parameter
#' @param auto Should the lengthscale and nugget values be automatically determined
#' by laGPs maximum likelihood algorithm? See \code{?laGP::mleGPsep} for more info
#' @param on_residuals Should a linear model be wrapped around the kriging model
#' to handle the main trends independently?
#'
#' @noRd
#' @keywords internal
interpolate <- function(independent, dependent, pred_grid, d = NA, g = NA, auto = F, on_residuals = T) {
  # check input
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  checkmate::assert_numeric(dependent, len = nrow(independent))
  checkmate::assert_class(pred_grid, "mobest_spatiotemporalpositions")
  if (on_residuals) {
    # linear fit
    combined <- independent %>% dplyr::mutate(d = dependent)
    model <- stats::lm(d ~ x + y + z, data = combined)
    dependent <- model[["residuals"]]
  }
  # priors for the global GP
  if ((any(is.na(d)) && is.na(g)) || auto) {
    da <- laGP::darg(list(mle = TRUE, max = 10), independent)
    ga <- laGP::garg(list(mle = TRUE, max = 10), dependent)
    d <- da$start
    g <- ga$start
  }
  # fit the global GP
  gp <- laGP::newGPsep(X = independent[, c("x", "y", "z")], Z = dependent, d = d, g = g)
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

#### helper functions ####

check_compatible_multi <- function(x, y, comp_f, ...) {
  purrr::pwalk(
    get_permutations(x, y),
    function(i1, i2) {
      comp_f(
        x[[i1]],
        y[[i2]],
        names(x)[i1],
        names(y)[i2],
        ...
      )
    }
  )
}

check_df_nrow_equal <- function(x, y, name_x, name_y) {
  if (nrow(x) != nrow(y)) {
    stop(name_x, " and ", name_y, ": Not the same number of lines")
  }
}

check_names_equal <- function(x, y, name_x, name_y, ignore_sd_cols = F) {
  names_x <- names(x)
  names_y <- names(y)
  if (ignore_sd_cols) {
    names_x <- names_x[!grepl("\\_sd$", names_x)]
    names_y <- names_y[!grepl("\\_sd$", names_y)]
  }
  if (!setequal(names_x, names_y)) {
    stop(name_x, " and ", name_y, ": Not the same names")
  }
}

get_permutations <- function(x, y, include.equals = FALSE) {
  expand.grid(seq_along(x), seq_along(y)) %>%
    as.data.frame() %>%
    magrittr::set_names(c("i1", "i2")) %>%
    dplyr::filter(
      .data[["i1"]] >= .data[["i2"]]
    )
}
