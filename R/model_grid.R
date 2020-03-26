#' create_model_grid
#'
#' @param independent_tables test
#' @param dependent_vars test
#' @param kernel_settings test
#'
#' @return test
#'
#' @export
create_model_grid <- function(independent_tables, dependent_vars, kernel_settings) {

  expand.grid(
    kernel_setting_id = kernel_settings$kernel_setting_id,
    dependent_var_id = c("PC1", "PC2", "PC3", "PC4"),
    independent_table_id = independent_tables$independent_table_id,
    stringsAsFactors = F
  ) %>%
    dplyr::left_join(
      kernel_settings, by = "kernel_setting_id"
    ) %>%
    dplyr::left_join(
      independent_tables, by = "independent_table_id"
    ) %>% dplyr::mutate(
      dependent_var = lapply(dependent_var_id, function(x) { anno[[x]] })
    ) %>% tibble::as_tibble()

}

#' run_model_grid
#'
#' @param model_grid test
#' @param pred_grid test
#'
#' @return test
#'
#' @export
run_model_grid <- function(model_grid, pred_grid) {

  # run interpolation for each entry in the model_grid
  prediction <- lapply(1:nrow(model_grid), function(i) {
    interpolate_laGP(
      independent = model_grid[["independent_table"]][[i]],
      dependent = model_grid[["dependent_var"]][[i]],
      pred_grid = pred_grid,
      auto = model_grid[["kernel_setting"]][[i]][["auto"]],
      d = model_grid[["kernel_setting"]][[i]][["d"]],
      g = model_grid[["kernel_setting"]][[i]][["g"]]
    )
  })

  # simplified model_grid
  model_grid_simplified <- model_grid %>%
    dplyr::mutate(independent_table_type = ifelse(independent_table_id == "age_center", "age_center", "age_sampled")) %>%
    dplyr::select(-kernel_setting, -independent_table, -dependent_var)

  # add prediction results for each run as a data.frame in a list column to model_grid
  model_grid_simplified$prediction_sample <- lapply(prediction, function(x) {
    pred <- data.frame(
      point_id = 1:length(x$mean),
      mean = x$mean,
      sd = sqrt(x$s2),
      stringsAsFactors = F
    )
    # merge with pred_grid to add relevant spatial information
    dplyr::left_join(
      pred,
      pred_grid,
      by = "point_id"
    )
  })

  return(model_grid_simplified)

}
