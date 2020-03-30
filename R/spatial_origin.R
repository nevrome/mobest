#' search_spatial_origin
#'
#' @param interpol_grid test
#'
#' @return test
#'
#' @export
search_spatial_origin <- function(interpol_grid, spatial_search_radius = 500000) {

  # transform runs for different PCs to columns
  pri <- tidyr::pivot_wider(
    interpol_grid,
    names_from = "dependent_var_id",
    values_from = c("mean", "sd")
  )

  # add new columns for output dataset
  pri <- dplyr::mutate(
      pri,
      angle = NA,
      genetic_distance = NA,
      spatial_distance = NA,
      mean_PC1_origin = NA,
      mean_PC2_origin = NA,
      mean_PC3_origin = NA,
      mean_PC4_origin = NA,
      x_origin = NA,
      y_origin = NA,
    )

  age_sample_run_pris <- split(pri, list(pri$independent_table_id, pri$kernel_setting_id))

  pri_ready_large <- pbapply::pblapply(age_sample_run_pris, function(age_sample_run_pri) {

    # split dataset by age slice
    time_pris <- age_sample_run_pri %>% split(age_sample_run_pri$z)

    for (p1 in 2:length(time_pris)) {

      # calculate spatial distance matrix between past and current points
      current_pri_spatial <- as.matrix(time_pris[[p1]][c("x", "y")])
      past_pri_spatial <- as.matrix(time_pris[[p1 - 1]][c("x", "y")])
      spatial_distance <- fields::rdist(current_pri_spatial, past_pri_spatial)

      # calculate PCA distance matrix between past and current points
      current_pri_genetics <- as.matrix(time_pris[[p1]][c("mean_PC1", "mean_PC2", "mean_PC3", "mean_PC4")])
      past_pri_genetics <- as.matrix(time_pris[[p1 - 1]][c("mean_PC1", "mean_PC2", "mean_PC3", "mean_PC4")])
      genetic_distance <- fields::rdist(current_pri_genetics, past_pri_genetics)

      # get points with least genetic distance in the past
      closest_point_indezes <- sapply(1:nrow(current_pri_genetics), function(x) {
        gendists <- genetic_distance[x,]
        gendists[spatial_distance[x,] > spatial_search_radius] <- NA
        which.min(gendists)
      })

      # add closest points info to current age slice points
      time_pris[[p1]] <- time_pris[[p1]] %>% dplyr::mutate(
        mean_PC1_origin = time_pris[[p1 - 1]]$mean_PC1[closest_point_indezes],
        mean_PC2_origin = time_pris[[p1 - 1]]$mean_PC2[closest_point_indezes],
        mean_PC3_origin = time_pris[[p1 - 1]]$mean_PC3[closest_point_indezes],
        mean_PC4_origin = time_pris[[p1 - 1]]$mean_PC4[closest_point_indezes],
        x_origin = time_pris[[p1 - 1]]$x[closest_point_indezes],
        y_origin = time_pris[[p1 - 1]]$y[closest_point_indezes]
      )
      time_pris[[p1]]$spatial_distance <- purrr::map2_dbl(
        1:length(closest_point_indezes), closest_point_indezes,
        function(i, j) { spatial_distance[i, j] }
      )
      time_pris[[p1]]$genetic_distance <- purrr::map2_dbl(
        1:length(closest_point_indezes), closest_point_indezes,
        function(i, j) { genetic_distance[i, j] }
      )

      # get spatial position of entangled points
      A <- as.matrix(time_pris[[p1]][c("x", "y")])
      B <- as.matrix(time_pris[[p1]][c("x_origin", "y_origin")])

      # calculate angle between points in radians and degrees
      AB <- B - A
      AC <- c(1, 0)
      time_pris[[p1]]$angle_rad <- sapply(
        1:nrow(time_pris[[p1]]), function(i) {
          if (time_pris[[p1]]$y_origin[i] < time_pris[[p1]]$y[i]) {
            2*pi - matlib::angle(AB[i,], AC, degree = FALSE)
          } else {
            matlib::angle(AB[i,], AC, degree = FALSE)
          }
        }
      )

      a_rad <- units::as_units(time_pris[[p1]]$angle_rad, "radians")
      a_deg <- units::set_units(a_rad, "degrees")
      time_pris[[p1]]$angle_degree <- as.numeric(a_deg)

    }

    # rowbind distance table
    pri_ready <- time_pris[2:length(time_pris)] %>% do.call(rbind, .)

    return(pri_ready)

  }, cl = 8)

  pri_ready <- pri_ready_large %>% dplyr::bind_rows()

  return(pri_ready)
}
