#' crossvalidate
#'
#' @param independent test
#' @param dependent test
#' @param kernel test
#' @param number_of_reorderings test
#'
#' @export
crossvalidate <- function(
  independent,
  dependent,
  kernel,
  number_of_reorderings = 1,
  number_of_splits = 10
) {

  crossval <- cbind(independent, dependent %>% dplyr::bind_cols())

  #### compile randomly reordered versions of crossval ####

  crossval_mixed_list <- lapply(1:number_of_reorderings, function(i) { crossval[sample(1:nrow(crossval), replace = F), ] })

  #### run prediction test for each of this versions ####

  lapply(crossval_mixed_list, function(crossval_mixed) {

    #### split crossval into 10 sections ####

    n <- number_of_splits
    nr <- nrow(crossval_mixed)
    crossval_10 <- split(crossval_mixed, rep(1:n, times = diff(floor(seq(0, nr, length.out = n + 1)))))

    # 9 sections are used as a training dataset for the GP model
    crossval_9_training <- lapply(
      1:n, function(i) {
        dplyr::bind_rows(crossval_10[-i])
      }
    )

    # 1 section is used as a test dataset
    crossval_9_test <- lapply(
      1:n, function(i) {
        crossval_10[[i]]
      }
    )

    #### prepare model grid for current (n-1):1 comparison with different kernels ####

    model_grid <- lapply(
      1:n, function(i) {

        # create model grid
        model_grid <- mobest::create_model_grid(
          independent = list(
            training = crossval_9_training[[i]] %>% dplyr::select(x, y, z)
          ),
          dependent = lapply(
            names(dependent), function(depvar) {
              crossval_9_training[[i]][[depvar]]
            }) %>%
            stats::setNames(names(dependent)),
          kernel = kernel,
          prediction_grid = list(
            crossval_9_test[[i]] %>%
              dplyr::select(x, y, z) %>%
              dplyr::mutate(point_id = 1:nrow(.))
          ) %>% stats::setNames(i)
        ) %>%
          dplyr::select(-independent_table_type)

      }
    ) %>% dplyr::bind_rows()

    #### run interpolation on model grid ####

    model_grid_result <- mobest::run_model_grid(model_grid)

    #### unnest prediction to get a point-wise prediction table ####

    interpol_grid <- mobest::unnest_model_grid(model_grid_result)

    #### merge prediction and real values ####

    # make wide for mean and sd PC values
    interpol_grid_wide <- interpol_grid %>% tidyr::pivot_wider(
      names_from = "dependent_var_id",
      values_from = c("mean", "sd")
    )

    # split by run training+test run setup to be able to merge with real test values
    interpol_grid_wide_split <- interpol_grid_wide %>% split(interpol_grid_wide$pred_grid_id)

    interpol_grid_merged <- lapply(
      as.numeric(unique(interpol_grid_wide$pred_grid_id)), function(i) {
        interpol_grid_wide_split[[i]] %>%
          dplyr::left_join(
            crossval_9_test[[i]] %>%
              dplyr::select(names(dependent)) %>%
              dplyr::mutate(point_id = 1:nrow(.)),
            by = "point_id"
          )
      }
    ) %>% dplyr::bind_rows()

    return(interpol_grid_merged)

  }) %>% dplyr::bind_rows() -> interpol_grid_merged_all

  return(interpol_grid_merged_all)

}

#' create_kernel_grid
#'
#' @param ds test
#' @param dt test
#' @param g test
#'
#' @export
create_kernel_grid <- function(ds, dt, g) {

  ks <- expand.grid(ds = ds, dt = dt, g = g)

  kernel_settings <- lapply(
      1:nrow(ks), function(i) {
        list(d = c(ks[["ds"]][i], ks[["ds"]][i], ks[["dt"]][i]), g = ks[["g"]][i], on_residuals = T, auto = F)
      }
    ) %>% setNames(
      sapply(
        1:nrow(ks), function(i) {
          paste0(ks[["ds"]][i]/1000, "_", ks[["dt"]][i], "_", ks[["g"]][i])
        }
      )
    )

  return(kernel_settings)

}
