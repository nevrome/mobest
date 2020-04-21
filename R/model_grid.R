#' create_model_grid
#'
#' @param data test
#' @param independent test
#' @param dependent test
#' @param kernel test
#' @param prediction_grid test
#'
#' @export
create_model_grid <- function(
  data,
  independent,
  dependent,
  kernel,
  prediction_grid
) {

  independent_tables <- tibble::tibble(
    independent_table = independent,
    independent_table_id = names(independent)
  )

  return(independent_tables)

  # # prep dependent vars
  # dependent_vars <- tibble::tibble(
  #   dependent_var_id = c("PC1", "PC2", "PC3", "PC4")
  # ) %>%
  #   dplyr::mutate(
  #     dependent_var = lapply(dependent_var_id, function(x) { anno[[x]] })
  #   )
  #
  # # create kernel parameters
  # kernel_settings <- tibble::tibble(
  #   kernel_setting = list(
  #     #ds50_dt100_g01 = list(auto = F, d = c(dist_scale_01_x_km(50), dist_scale_01_x_km(50), dist_scale_01_z_years(100)), g = 0.1),
  #     #ds100_dt200_g01 = list(auto = F, d = c(dist_scale_01_x_km(100), dist_scale_01_x_km(100), dist_scale_01_z_years(200)), g = 0.1),
  #     #ds500_dt500_g01 = list(d = c(500000, 500000, 500), g = 0.1, on_residuals = T, auto = F),
  #     ds1000_dt1000_g01 = list(d = c(1000000, 1000000, 1000), g = 0.1, on_residuals = T, auto = F)
  #   ),
  #   kernel_setting_id = names(kernel_setting)
  # )
  #
  # # individual point
  # sf::st_as_sf(
  #   tibble::tibble(lon = 9.05, lat = 48.52),
  #   coords = c("lon", "lat"),
  #   crs = 4326,
  #   remove = FALSE
  # ) %>%
  #   sf::st_transform(
  #     crs = "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs",
  #   ) %>% sf::st_coordinates()
  #
  # # create spatiotemporal prediction grid
  # pred_grids <- tibble::tibble(
  #   pred_grid = list(
  #     #scs100_tl100 = mobest::create_prediction_grid(area, spatial_cell_size = 100000, time_layers = seq(-7500, -500, 100)),
  #     #scs200_tl200 = mobest::create_prediction_grid(area, spatial_cell_size = 200000, time_layers = seq(-7500, -500, 200))
  #     tuebingen = tibble::tibble(x = -69459.46, y = 2031623, z = seq(-7500, -500, 100), point_id = 1:71)
  #   ),
  #   pred_grid_id = names(pred_grid)
  # )
  #
  # # merge info in prepare model grid
  # model_grid_pca <- mobest::create_model_grid(
  #   independent_tables = independent_tables,
  #   dependent_vars = dependent_vars,
  #   kernel_settings = kernel_settings,
  #   pred_grids = pred_grids
  # )
  #
  # #### prepare mds model grid ####
  #
  # anno_mds <- anno %>% dplyr::filter(
  #   !is.na(C1)
  # )
  #
  # # prep independent variables with temporal sampling
  # independent_tables <- tibble::tibble(
  #   independent_table = c(
  #     list(dplyr::transmute(.data = anno_mds, x = x, y = y, z = calage_center))
  #   ),
  #   independent_table_id = c("age_center")
  # )
  #
  # # prep dependent vars
  # dependent_vars <- tibble::tibble(
  #   dependent_var_id = c("C1", "C2", "C3", "C4")
  # ) %>%
  #   dplyr::mutate(
  #     dependent_var = lapply(dependent_var_id, function(x) { anno_mds[[x]] })
  #   )
  #
  # # merge info in prepare model grid
  # model_grid_mds <- mobest::create_model_grid(
  #   independent_tables = independent_tables,
  #   dependent_vars = dependent_vars,
  #   kernel_settings = kernel_settings,
  #   pred_grids = pred_grids
  # )

}


#' create_model_grid_raw
#'
#' @param independent_tables test
#' @param dependent_vars test
#' @param kernel_settings test
#' @param pred_grids test
#'
#' @return test
#'
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
    ) %>%
    dplyr::mutate(
      independent_table_type = ifelse(.data[["independent_table_id"]] == "age_center", "age_center", "age_sampled")
    ) %>%
    tibble::as_tibble()

}

#' run_model_grid
#'
#' @param model_grid test
#' @param quiet test
#'
#' @return test
#'
#' @export
run_model_grid <- function(model_grid, quiet = F) {

  # run interpolation for each entry in the model_grid
  prediction <- lapply(1:nrow(model_grid), function(i) {
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
  })

  # simplify model_grid
  model_grid_simplified <- model_grid %>%
    dplyr::select(
      -.data[["kernel_setting"]],
      -.data[["independent_table"]],
      -.data[["dependent_var"]],
      -.data[["pred_grid"]]
    )

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
