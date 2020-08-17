#' search_spatial_origin
#'
#' @param interpol_grid test
#' @param spatial_search_radius test
#'
#' @return test
#'
#' @export
search_spatial_origin <- function(interpol_grid, spatial_search_radius = 500000) {

  # transform runs for different ancestry components to columns
  pri <- tidyr::pivot_wider(
    interpol_grid,
    names_from = "dependent_var_id",
    values_from = c("mean", "sd")
  )

  # add new columns for output dataset
  pri <- dplyr::mutate(
      pri,
      angle_deg = NA,
      genetic_distance = NA,
      spatial_distance = NA,
      mean_C1_origin = NA,
      mean_C2_origin = NA,
      x_origin = NA,
      y_origin = NA,
    )

  age_sample_run_pris <- split(pri, list(pri[["independent_table_id"]], pri[["kernel_setting_id"]], pri[["pred_grid_id"]]))

  pri_ready_large <- pbapply::pblapply(age_sample_run_pris, function(age_sample_run_pri) {

    # split dataset by age slice
    time_pris <- age_sample_run_pri %>% split(age_sample_run_pri[["z"]])

    for (p1 in 2:length(time_pris)) {

      # calculate spatial distance matrix between past and current points
      current_pri_spatial <- as.matrix(time_pris[[p1]][c("x", "y")])
      past_pri_spatial <- as.matrix(time_pris[[p1 - 1]][c("x", "y")])
      spatial_distance <- fields::rdist(current_pri_spatial, past_pri_spatial)

      # calculate genetic distance matrix between past and current points
      current_pri_genetics <- as.matrix(time_pris[[p1]][c("mean_C1", "mean_C2")])
      past_pri_genetics <- as.matrix(time_pris[[p1 - 1]][c("mean_C1", "mean_C2")])
      genetic_distance <- fields::rdist(current_pri_genetics, past_pri_genetics)

      # get points with least genetic distance in the past
      closest_point_indezes <- sapply(1:nrow(current_pri_genetics), function(x) {
        # all genetic distances to current point
        gendists <- genetic_distance[x,]
        # reduce selection of genetic distances to only include distances within the spatial distance radius
        gendists[spatial_distance[x,] > spatial_search_radius] <- NA
        # find all points with min genetic distances
        min_gen_distance_points <- which(gendists == min(gendists, na.rm = TRUE))
        # when multiple: return the spatially closest. If there are multiple equally close
        # in space then this will return the first of these points, as which.min always
        # just returns the first element in case of multiple solutions
        if (length(min_gen_distance_points) > 1) {
          min_gen_distance_points[which.min(spatial_distance[x, ][min_gen_distance_points])]
        } else {
          min_gen_distance_points
        }
      })

      # add closest points info to current age slice points
      time_pris[[p1]] <- time_pris[[p1]] %>% dplyr::mutate(
        mean_C1_origin = time_pris[[p1 - 1]]$mean_C1[closest_point_indezes],
        mean_C2_origin = time_pris[[p1 - 1]]$mean_C2[closest_point_indezes],
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

    }

    # rowbind distance table
    pri_ready <- time_pris[2:length(time_pris)] %>% do.call(rbind, .)

    return(pri_ready)

  }, cl = parallel::detectCores())

  pri_ready <- pri_ready_large %>% dplyr::bind_rows()

  # add add_origin_vector_coordinates
  pri_ready <- pri_ready %>% add_origin_vector_coordinates()

  # add angle
  pri_ready$angle_deg[pri_ready$spatial_distance != 0] <- sapply(
    1:nrow(pri_ready[pri_ready$spatial_distance != 0, ]), function(i) {
    vec2deg(c(pri_ready$x_to_origin[i], pri_ready$y_to_origin[i]))
  })

  return(pri_ready)
}

add_origin_vector_coordinates <- function(x) {

  x <- x %>%
    dplyr::mutate(
      x_to_origin = .data[["x_origin"]] - .data[["x"]],
      y_to_origin = .data[["y_origin"]] - .data[["y"]]
    )

  normalized_vector <- purrr::map2(
    x[["x_to_origin"]], x[["y_to_origin"]], function(l, r) {
      scalar1(c(l, r))
    }
  ) %>% do.call(rbind, .)

  x <- x %>%
    dplyr::mutate(
      x_to_origin_norm = normalized_vector[,1],
      y_to_origin_norm = normalized_vector[,2]
    )

  return(x)
}

scalar1 <- function(x) {
  if (all(x == 0)) {
    x
  } else {
    x / sqrt(sum(x^2))
  }
}

