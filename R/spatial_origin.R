#' search_spatial_origin
#'
#' @param interpol_grid test
#' @param nugget test
#'
#' @return test
#'
#' @export
search_spatial_origin <- function(interpol_grid, nugget = 0.01) {

  dependent_vars <- unique(interpol_grid$dependent_var_id)
  mean_cols <- paste0("mean_", dependent_vars)
  sd_cols <- paste0("sd_", dependent_vars)

  # remove prediction points with too high standard deviation
  interpol_grid_sd_filtered <- interpol_grid %>%
    dplyr::group_by(independent_table_id, dependent_var_id, kernel_setting_id, pred_grid_id) %>%
    dplyr::filter(
      sd < 0.2 * diff(range(mean))
    ) %>%
    dplyr::ungroup()

  # transform runs for different ancestry components to columns
  pri <- tidyr::pivot_wider(
    interpol_grid_sd_filtered,
    names_from = "dependent_var_id",
    values_from = c("mean", "sd")
  )

  # filter prediction points that were only half removed by the sd filter
  pri <- pri %>%
    dplyr::filter(
      dplyr::across(tidyselect::one_of(mean_cols), ~!is.na(.x))
    )

  # add new columns for output dataset
  pri <- pri %>% dplyr::mutate(
      spatial_distance = NA_real_,
      x_origin = NA_real_,
      y_origin = NA_real_,
      angle_deg = NA_real_,
    )

  # split by
  age_sample_run_pris <- pri %>% dplyr::group_split(
    independent_table_id, kernel_setting_id, pred_grid_id
  )

  # loop by
  pri_ready_large <- pbapply::pblapply(age_sample_run_pris, function(age_sample_run_pri) {

    # split dataset by age slice
    time_pris <- split(
      age_sample_run_pri,
      age_sample_run_pri[["z"]]
    )

    for (p1 in 2:length(time_pris)) {

      # calculate spatial distance matrix between past and current points
      current_pri_spatial <- as.matrix(time_pris[[p1]][c("x", "y")])
      past_pri_spatial <- as.matrix(time_pris[[p1 - 1]][c("x", "y")])
      spatial_distance <- fields::rdist(current_pri_spatial, past_pri_spatial)

      # calculate genetic distance matrix between past and current points
      current_pri_genetics <- as.matrix(time_pris[[p1]][mean_cols])
      past_pri_genetics <- as.matrix(time_pris[[p1 - 1]][mean_cols])
      past_pri_genetics_sd <- as.matrix(time_pris[[p1 - 1]][sd_cols])
      genetic_distance <- fields::rdist(current_pri_genetics, past_pri_genetics)

      # get points with least genetic distance in the past
      centroid_points <- do.call(rbind, lapply(1:nrow(current_pri_genetics), function(index_of_A) {
        # all genetic distances to current point A
        gendists_to_A <- genetic_distance[index_of_A,]
        # find closest point in the past B
        index_of_B <- which.min(gendists_to_A)
        # find points with similar genetic makeup like B
        B_mean <- past_pri_genetics[index_of_B,]
        B_sd <- past_pri_genetics_sd[index_of_B,]
        B_spatial_points <- past_pri_spatial[
          past_pri_genetics[,1] < B_mean[1] + nugget & past_pri_genetics[,1] > B_mean[1] - nugget &
          past_pri_genetics[,2] < B_mean[2] + nugget & past_pri_genetics[,2] > B_mean[2] - nugget,
        ]
        # find centroid point C
        if (is.vector(B_spatial_points)) {
          C <- c(B_spatial_points[1], B_spatial_points[2])
        } else {
          C <- c(mean(B_spatial_points[,1]), mean(B_spatial_points[,2]))
        }
        # return centroid point
        return(C)
      }))

      # add closest points info to current age slice points
      time_pris[[p1]]$x_origin <- centroid_points[,1]
      time_pris[[p1]]$y_origin <- centroid_points[,2]
    }

    # rowbind distance table
    pri_ready <- time_pris[2:length(time_pris)] %>% do.call(rbind, .)

    return(pri_ready)

  }, cl = parallel::detectCores())

  pri_ready <- pri_ready_large %>% dplyr::bind_rows()

  # add add_origin_vector_coordinates
  pri_ready <- pri_ready %>% add_origin_vector_coordinates()

  # add distance
  pri_ready$spatial_distance <- sqrt(pri_ready$x_to_origin^2 + pri_ready$y_to_origin^2)

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

