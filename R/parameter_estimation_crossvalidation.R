#' crossvalidate
#'
#' @param independent test
#' @param dependent test
#' @param kernel test
#' @param iterations test
#' @param groups test
#' @param quiet
#'
#' @export
crossvalidate <- function(
  independent,
  dependent,
  kernel,
  iterations,
  groups = 10,
  quiet = T
) {
  # input check
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  checkmate::assert_class(dependent, "mobest_observations")
  checkmate::assert_list(
    kernel, types = "mobest_kernelsetting",
    any.missing = F, min.len = 1, names = "strict"
  )
  checkmate::assert_count(iterations)
  checkmate::assert_count(groups)
  # create crossvalidation dataset
  crossval <- cbind(independent, dependent %>% dplyr::bind_cols())
  # compile randomly reordered versions of crossval
  crossval_mixed_list <- lapply(1:iterations, function(i) {
    dplyr::slice_sample(crossval, n = nrow(crossval), replace = F)
  })
  # run prediction test for each iteration
  crossval_interpol_grid <- purrr::map2_dfr(
    1:iterations, crossval_mixed_list,
    function(mixing_iteration, crossval_mixed) {
      # split crossval into sections
      n <- groups
      nr <- nrow(crossval_mixed)
      crossval_all <- split(crossval_mixed, rep(1:n, times = diff(floor(seq(0, nr, length.out = n + 1)))))
      # n-1 sections are used as a training dataset for the GP model
      crossval_training <- purrr::map(1:n, function(i) { dplyr::bind_rows(crossval_all[-i]) })
      # 1 section is used as a test dataset
      crossval_test <- purrr::map(1:n, function(i) { crossval_all[[i]] })
      # prepare model grid for current (n-1):1 comparison with different kernels
      model_grid <- suppressMessages(
        purrr::pmap_dfr(
          list(1:n, crossval_training, crossval_test),
          function(run_id, training, test) {
            mobest::create_model_grid(
              independent = mobest::create_spatpos_multi(
                id = training$id,
                x = list(training$x),
                y = list(training$y),
                z = list(training$z),
                it = paste0("ind_crossval_run_", run_id)
              ),
              dependent = do.call(
                mobest::create_obs,
                training[, names(dependent)]
              ),
              kernel = kernel,
              prediction_grid = mobest::create_spatpos_multi(
                id = test$id,
                x = list(test$x),
                y = list(test$y),
                z = list(test$z),
                it = paste0("pred_crossval_run_", run_id)
              )
            )
          }
        )
      )
      # run interpolation on model grid
      interpol_grid <- mobest::run_model_grid(model_grid, quiet = quiet)
      # add mixing iteration column
      interpol_grid %>% dplyr::mutate(mixing_iteration = mixing_iteration)
    }
  )
  # make wide for predicted (by the GPR) mean and sd values
  crossval_interpol_grid_wide <- crossval_interpol_grid %>% tidyr::pivot_wider(
    names_from = "dependent_var_id",
    values_from = c("mean", "sd")
  )
  # merge with actually measured information
  crossval_interpol_comparison <- crossval_interpol_grid_wide %>%
    dplyr::left_join(
      crossval %>% dplyr::select(
        -.data[["x"]], -.data[["y"]], -.data[["z"]]
      ),
      by = "id"
    )
  # calculate differences between estimated and measured values
  for (dep in names(dependent)) {
    crossval_interpol_comparison[[paste0(dep, "_dist")]] <-
      crossval_interpol_comparison[[paste0("mean_", dep)]] -
      crossval_interpol_comparison[[dep]]
  }
  # prepare output dataset
  crossval_interpol_comparison %>%
    dplyr::select(
      .data[["id"]],
      .data[["mixing_iteration"]],
      .data[["kernel_setting_id"]],
      tidyselect::contains("_dist")
    ) %>%
    tidyr::pivot_longer(
      cols = tidyselect::contains("_dist"),
      names_to = "dependent_var",
      values_to = "difference"
    ) %>%
    # turn kernel parameters into distinct columns again
    dplyr::mutate(
      kernel_setting_id = gsub("kernel_", "", .data[["kernel_setting_id"]])
    ) %>%
    tidyr::separate(
      .data[["kernel_setting_id"]],
      c("ds", "dt", "g"),
      sep = "_",
      convert = T,
      remove = F
    ) %>%
    dplyr::select(
      -.data[["kernel_setting_id"]]
    )
}
