#' crossvalidate
#'
#' @param independent test
#' @param dependent test
#' @param kernel test
#' @param iterations test
#' @param groups test
#'
#' @export
crossvalidate <- function(
  independent,
  dependent,
  kernel,
  iterations,
  groups = 10
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
  purrr::map_dfr(crossval_mixed_list, function(crossval_mixed) {
    # split crossval into sections
    n <- groups
    nr <- nrow(crossval_mixed)
    crossval_all <- split(crossval_mixed, rep(1:n, times = diff(floor(seq(0, nr, length.out = n + 1)))))
    # n-1 sections are used as a training dataset for the GP model
    crossval_training <- purrr::map(1:n, function(i) { dplyr::bind_rows(crossval_all[-i]) })
    # 1 section is used as a test dataset
    crossval_test <- purrr::map(1:n, function(i) { crossval_10[[i]] })
    # prepare model grid for current (n-1):1 comparison with different kernels
    model_grid <- purrr::map2_dfr(
      crossval_training, crossval_test,
      function(training, test) {
        mobest::create_model_grid(
          independent = mobest::create_spatpos_multi(
            id = training$id,
            x = list(training$x),
            y = list(training$y),
            z = list(training$z),
            it = "age_median"
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
            it = "age_median"
          )
        )
      }
    )
    # run interpolation on model grid
    interpol_grid <- mobest::run_model_grid(model_grid)
    # merge prediction and real values
    # make wide for mean and sd PC values
    interpol_grid_wide <- interpol_grid %>% tidyr::pivot_wider(
      names_from = "dependent_var_id",
      values_from = c("mean", "sd")
    )

    #### TODO from here ####

    # split by run training+test run setup to be able to merge with real test values
    interpol_grid_wide_split <- interpol_grid_wide %>% split(interpol_grid_wide$pred_grid_id)

    interpol_grid_merged <- lapply(
      unique(interpol_grid_wide$pred_grid_id), function(i) {
        interpol_grid_wide_split[[i]] %>%
          dplyr::left_join(
            crossval_9_test[[as.numeric(i)]] %>%
              dplyr::select(names(dependent)) %>%
              dplyr::mutate(point_id = 1:nrow(.)),
            by = "point_id"
          )
      }
    ) %>% dplyr::bind_rows()

    return(interpol_grid_merged)

  }) -> interpol_grid_merged_all

  for (dep in names(dependent)) {
    interpol_grid_merged_all[[paste0(dep, "_dist")]] <- interpol_grid_merged_all[[dep]] -
      interpol_grid_merged_all[[paste0("mean_", dep)]]
  }

  interpol_comparison <- interpol_grid_merged_all %>%
    dplyr::select(
      kernel_setting_id, tidyselect::contains("_dist")
    ) %>%
    tidyr::pivot_longer(
      cols = tidyselect::contains("_dist"),
      names_to = "dependent_var",
      values_to = "difference"
    )

  # turn kernel parameters into distinct columns again
  interpol_comparison <- interpol_comparison %>%
    tidyr::separate(
      kernel_setting_id,
      c("ds", "dt", "g"),
      sep = "_",
      convert = T,
      remove = F
    )

  return(interpol_comparison)

}
