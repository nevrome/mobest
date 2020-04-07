#' create_model_grid
#'
#' @param independent_tables test
#' @param dependent_vars test
#' @param kernel_settings test
#' @param pred_grids test
#'
#' @return test
#'
#' @export
create_model_grid <- function(independent_tables, dependent_vars, kernel_settings, pred_grids) {

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
    ) %>%
    dplyr::mutate(
      independent_table_type = ifelse(independent_table_id == "age_center", "age_center", "age_sampled")
    ) %>%
    tibble::as_tibble()

}

#' run_model_grid
#'
#' @param model_grid test
#'
#' @return test
#'
#' @export
run_model_grid <- function(model_grid) {

  # run interpolation for each entry in the model_grid
  prediction <- lapply(1:nrow(model_grid), function(i) {
    interpolate_laGP(
      independent = model_grid[["independent_table"]][[i]],
      dependent = model_grid[["dependent_var"]][[i]],
      pred_grid = model_grid[["pred_grid"]][[i]],
      auto = model_grid[["kernel_setting"]][[i]][["auto"]],
      d = model_grid[["kernel_setting"]][[i]][["d"]],
      g = model_grid[["kernel_setting"]][[i]][["g"]]
    )
  })

  # simplify model_grid
  model_grid_simplified <- model_grid %>%
    dplyr::select(-kernel_setting, -independent_table, -dependent_var, -pred_grid)

  # add prediction results for each run as a data.frame in a list column to model_grid
  model_grid_simplified$prediction_sample <- purrr::map2(prediction, model_grid[["pred_grid"]], function(x, y) {
    pred <- data.frame(
      point_id = 1:length(x$mean),
      mean = x$mean,
      sd = sqrt(x$s2),
      stringsAsFactors = F
    )
    # merge with pred_grid to add relevant spatial information
    dplyr::left_join(
      pred, y,
      by = "point_id"
    )
  })

  return(model_grid_simplified)
}
