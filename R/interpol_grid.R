#' unnest_model_grid
#'
#' @param model_grid test
#'
#' @return test
#'
#' @export
unnest_model_grid <- function(model_grid) {

  model_grid %>%
    tidyr::unnest(cols = "prediction_sample")

}

#' condense_interpol_grid
#'
#' @param interpol_grid test
#'
#' @return test
#'
#' @export
condense_interpol_grid <- function(interpol_grid) {

  interpol_grid %>%
    dplyr::group_by(x, y, z, point_id, independent_table_type, dependent_var_id, kernel_setting_id, pred_grid_id) %>%
    dplyr::summarize(
      sd = age_center_catering_sd(independent_table_type, mean, sd),
      mean = mean(mean)
    ) %>%
    dplyr::ungroup()

}

age_center_catering_sd <- function(independent_table_type, input_mean, input_sd) {
  if (unique(independent_table_type) == "age_center") {
    input_sd
  } else {
    sd(sapply(1:length(input_mean), function(i) { rnorm(1, input_mean[i], input_sd[i]) }))
  }
}
