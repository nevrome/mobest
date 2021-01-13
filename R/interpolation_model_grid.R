#' Create a kriging model grid
#'
#' Constructs a model grid with all combinations of the different input parameter
#' configurations for the kriging model.
#'
#' @param independent Named list of dataframes. Spatiotemporal input point positions.
#' Each dataframe should have three numeric columns x, y and z:
#'
#' \itemize{
#'  \item{x: }{Spatial coordinate in x-axis direction (in a cartesian grid)}
#'  \item{y: }{Spatial coordinate in y-axis direction}
#'  \item{z: }{Temporal position (age)}
#' }
#'
#' @param dependent Named list of numeric vectors.
#' Dependent variables that should be interpolated. Each vector should have one
#' entry for each row in the \code{independent} list dataframes
#' @param kernel Named list of lists. Kernel parameter settings:
#'
#' \itemize{
#'  \item{: }{Numeric vector with lengthscale values}
#'  \item{g: }{Nugget value}
#'  \item{on_residuals: }{Should a linear model take out the main trends before the kriging interpolation?}
#'  \item{auto: }{Should the lengthscale and nugget values be automatically determined by laGPs maximum likelihood algorithm?}
#' }
#'
#' See \code{?interpolate_laGP} for more information
#'
#' @param prediction_grid Named list of dataframes.
#' Prediction grid positions for the interpolation.
#' Each dataframe should have three numeric columns x, y and z and a point id column:
#'
#' \itemize{
#'  \item{point_id: }{Unique point id}
#'  \item{x: }{Spatial coordinate in x-axis direction (in a cartesian grid)}
#'  \item{y: }{Spatial coordinate in y-axis direction}
#'  \item{z: }{Temporal position (age)}
#'  \item{region_id: }{Name of the region this point belongs to}
#' }
#'
#' See \code{?create_prediction_grid} for a function to create grid for a certain
#' spatial region
#'
#' @return An object of class \code{mobest_modelgrid} which inherits from tibble
#'
#' @export
create_model_grid <- function(
  independent,
  dependent,
  kernel,
  prediction_grid
) {
  # input check
  checkmate::assert_list(
    independent, types = "mobest_spatiotemporalpositions",
    any.missing = F, min.len = 1, names = "strict"
  )
  checkmate::assert_list(
    dependent, types = "numeric",
    any.missing = F, min.len = 1, names = "strict"
  )
  checkmate::assert_list(
    kernel, types = "mobest_kernelsetting",
    any.missing = F, min.len = 1, names = "strict"
  )
  checkmate::assert_list(
    prediction_grid, types = "mobest_spatiotemporalpositions",
    any.missing = F, min.len = 1, names = "strict"
  )
  # fill create general structure and id columns
  independent_tables <- tibble::tibble(
    independent_table = independent,
    independent_table_id = factor(names(independent), levels = names(independent))
  )
  dependent_vars <- tibble::tibble(
    dependent_var = dependent,
    dependent_var_id = factor(names(dependent), levels = names(dependent))
  )
  kernel_settings <- tibble::tibble(
    kernel_setting = kernel,
    kernel_setting_id = factor(names(kernel), levels = names(kernel))
  )
  pred_grids <- tibble::tibble(
    pred_grid = prediction_grid,
    pred_grid_id = factor(names(prediction_grid), levels = names(prediction_grid))
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
    independent_table_id = independent_tables[["independent_table_id"]],
    dependent_var_id = dependent_vars[["dependent_var_id"]],
    kernel_setting_id = kernel_settings[["kernel_setting_id"]],
    pred_grid_id = pred_grids[["pred_grid_id"]],
    stringsAsFactors = F
  ) %>%
    dplyr::left_join(
      independent_tables, by = "independent_table_id"
    ) %>%
    dplyr::left_join(
      dependent_vars, by = "dependent_var_id"
    ) %>%
    dplyr::left_join(
      kernel_settings, by = "kernel_setting_id"
    ) %>%
    dplyr::left_join(
      pred_grids, by = "pred_grid_id"
    )
}

#' Run kriging interpolation on model grid
#'
#' Calculate kriging interpolation for all entries in a model grid.
#'
#' @param model_grid An object of class \code{mobest_modelgrid} as created by
#' \link{create_model_grid}
#' @param unnest Boolean. Should the kriging result be unnested to return a
#' prediction point-wise table of class \code{mobest_interpolgrid}?
#' @param quiet Boolean. Should a progress indication be printed?
#'
#' @return If \code{unnest = T } then an object of class \code{mobest_interpolgrid},
#' otherwise a tibble with a list column \code{prediction} that contains the
#' kriging results for each model grid row
#'
#' @rdname run_model_grid
#' @export
run_model_grid <- function(model_grid, unnest = T, quiet = F) {
  UseMethod("run_model_grid")
}

#' @rdname run_model_grid
#' @export
run_model_grid.default <- function(model_grid, unnest = T, quiet = F) {
  stop("x is not an object of class mobest_modelgrid")
}

#' @rdname run_model_grid
#' @export
run_model_grid.mobest_modelgrid <- function(model_grid, unnest = T, quiet = F) {
  # run interpolation for each entry in the model_grid
  prediction <- purrr::map(
    1:nrow(model_grid), function(i) {
      if (!quiet) {
        message("running model ", i, " of ", nrow(model_grid))
      }
      interpolate_laGP(
        independent = model_grid[["independent_table"]][[i]],
        dependent = model_grid[["dependent_var"]][[i]],
        pred_grid = model_grid[["pred_grid"]][[i]],
        # d has to be squared because of the configuration of the default laGP kernel
        d = model_grid[["kernel_setting"]][[i]][["d"]]^2,
        g = model_grid[["kernel_setting"]][[i]][["g"]],
        auto = model_grid[["kernel_setting"]][[i]][["auto"]],
        on_residuals = model_grid[["kernel_setting"]][[i]][["on_residuals"]]
      )
    }
  )
  # simplify model_grid
  model_grid_simplified <- model_grid %>%
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
