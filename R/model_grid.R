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
