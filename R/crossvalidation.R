#' Use crossvalidation for kriging kernel parameter estimation
#'
#' @param independent An object of class mobest_spatiotemporalpositions
#' @param dependent An object of class mobest_observations
#' @param kernel An object of class mobest_kernelsetting
#' @param iterations Integer. Number of crossvalidation iterations. Each iteration
#' goes along with a random reordering of the input points for training and test
#' data
#' @param groups Integer. Number of groups for splitting up training and test data
#' 10 means for example that the data should be split up into 9 parts training and
#' 1 part test observations
#' @param quiet Logical. Should a progress indication be printed?
#'
#' @export
crossvalidate <- function(
  independent,
  dependent,
  kernel,
  iterations,
  groups = 10,
  quiet = F
) {
  # input check
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  checkmate::assert_class(dependent, "mobest_observations")
  checkmate::assert_class(kernel, "mobest_kernelsetting_multi")
  checkmate::assert_count(iterations)
  checkmate::assert_count(groups)
  # create crossvalidation dataset
  if (!quiet) { message("Preparing crossvalidation dataset") }
  crossval <- cbind(independent, dependent %>% dplyr::bind_cols())
  # compile randomly reordered versions of crossval
  crossval_mixed_list <- lapply(1:iterations, function(i) {
    dplyr::slice_sample(crossval, n = nrow(crossval), replace = F)
  })
  # run prediction test for each iteration
  if (!quiet) { message("Running mixing iterations") }
  crossval_interpol_grid <- purrr::map2_dfr(
    1:iterations, # this counter is only passed here to document the run number in the output df
    crossval_mixed_list,
    function(mixing_iteration, crossval_mixed) {
      if (!quiet) { message("Starting mixing iteration ", mixing_iteration, " of ", iterations) }
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
                mobest::create_spatpos(
                  id = training$id,
                  x = training$x, y = training$y, z = training$z
                ),
                .names = paste0("ind_crossval_run_", run_id)
              ),
              dependent = mobest::create_obs_multi(
                do.call(
                  mobest::create_obs,
                  training[names(dependent)]
                ),
                .names = paste0("obs_crossval_run_", run_id)
              ),
              kernel = kernel,
              prediction_grid = mobest::create_spatpos_multi(
                mobest::create_spatpos(
                  id = test$id,
                  x = test$x, y = test$y, z = test$z
                ),
                .names = paste0("pred_crossval_run_", run_id)
              )
            )
          }
        )
      )
      # run interpolation on model grid
      interpol_grid <- mobest::run_model_grid(model_grid, quiet = quiet)
      interpol_grid %>%
        # add mixing iteration column
        dplyr::mutate(
          mixing_iteration = mixing_iteration,
          .after = "pred_grid_id"
        )
    }
  )
  # merge with actually measured information
  if (!quiet) { message("Comparing estimated with measured information") }
  crossval_interpol_grid %>%
    dplyr::left_join(
      crossval %>% dplyr::select(
        -.data[["x"]], -.data[["y"]], -.data[["z"]]
      ) %>% tidyr::pivot_longer(
        cols = -"id",
        names_to = "dependent_var_id",
        values_to = "measured"
      ),
      by = c("id", "dependent_var_id")
    ) %>%
    # calculate differences between estimated and measured values
    dplyr::mutate(
      difference = .data[["mean"]] - .data[["measured"]]
    )
}
