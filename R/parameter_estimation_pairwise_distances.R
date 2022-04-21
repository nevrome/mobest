#' Calculate pairwise distances and prepare an empirical semivariogram
#'
#' @param independent An object of class mobest_spatiotemporalpositions
#' @param dependent An object of class mobest_observations
#' @param m_to_km Logical. Should distances be transformed from m to km (x/1000)
#' @param with_resid Logical. Calculate distances also on the residuals of a spatiotemporal linear model
#' @param pairwise_distances An object of class mobest_pairwisedistances
#' @param geo_bin Numeric. Width of the spatial bins
#' @param time_bin Numeric. Width of the temporal bins
#' @param per_bin_operation Function. Summarising operation that should be applied to each distance.
#' Default: output = half mean squared input.
#'
#' @rdname pairwise_distances
#' @export
calculate_pairwise_distances <- function(independent, dependent, m_to_km = T, with_resid = T) {
  # input check
  if (nrow(independent) != nrow(dependent)) {
    stop("independent and dependent must have the same number of rows")
  }
  # geo distance
  d_geo_long <- calculate_geo_pairwise_distances(independent, m_to_km = m_to_km)
  # time distance
  d_time_long <- calculate_time_pairwise_distances(independent)
  # obs total distance
  d_obs_total <- stats::dist(dependent, "euclidean") %>% as.matrix()
  rownames(d_obs_total) <- colnames(d_obs_total) <- independent[["id"]]
  d_obs_total_long <- d_obs_total %>%
    reshape2::melt(value.name = "obs_dist_total") %>%
    dplyr::rename(id1 = "Var1", id2 = "Var2")
  # obs individual distance
  d_obs_long_list <- calculate_dependent_pairwise_distances(
    independent[["id"]], dependent, with_resid = with_resid, independent
  )
  # join different distances
  purrr::reduce(
    list(d_geo_long, d_time_long, d_obs_total_long, d_obs_long_list),
    function(x, y) {
      dplyr::full_join(
        x, y, by = c("id1", "id2")
      )
    }
  ) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_pairwisedistances")
}

#' @rdname pairwise_distances
#' @export
calculate_geo_pairwise_distances <- function(independent, m_to_km = T) {
  # input checks
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  # calculate distances
  d_geo <- stats::dist(independent %>% dplyr::select(.data[["x"]], .data[["y"]]), "euclidean") %>% as.matrix()
  rownames(d_geo) <- colnames(d_geo) <- independent[["id"]]
  d_geo %>%
    reshape2::melt(value.name = "geo_dist") %>%
    dplyr::rename(id1 = "Var1", id2 = "Var2") %>%
    dplyr::mutate(
      # m to km
      geo_dist = if (m_to_km) {.data[["geo_dist"]]/1000} else {.data[["geo_dist"]]}
    ) %>%
    tibble::as_tibble()
}

#' @rdname pairwise_distances
#' @export
calculate_time_pairwise_distances <- function(independent) {
  # input checks
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  # calculate distances
  d_time <- stats::dist(independent %>% dplyr::select(.data[["z"]]), "euclidean") %>% as.matrix()
  rownames(d_time) <- colnames(d_time) <- independent[["id"]]
  d_time %>%
    reshape2::melt(value.name = "time_dist") %>%
    dplyr::rename(id1 = "Var1", id2 = "Var2") %>%
    tibble::as_tibble()
}

#' @rdname pairwise_distances
#' @export
calculate_dependent_pairwise_distances <- function(ids, dependent, with_resid = F, independent = NULL) {
  # input checks
  checkmate::assert_class(dependent, "mobest_observations")
  if (with_resid & is.null(independent)) { stop("If with_resid, then independent can not be NULL") }
  # calculate distances
  var_names <- names(dependent)
  d_obs_long_list <- var_names %>%
    purrr::map(
      function(var_name) {
        # d_obs
        d_obs <- stats::dist(dependent[[var_name]], "euclidean") %>% as.matrix()
        rownames(d_obs) <- colnames(d_obs) <- ids
        d_obs_long <- d_obs %>%
          reshape2::melt(value.name = paste0(var_name, "_dist")) %>%
          dplyr::rename(id1 = "Var1", id2 = "Var2")
        # d_obs_resid
        if (with_resid) {
          model <- stats::lm(dependent[[var_name]] ~ independent$x + independent$y + independent$z)
          model_residuals <- stats::residuals(model)
          d_obs_resid <- as.matrix(stats::dist(model_residuals, "euclidean"))
          rownames(d_obs_resid) <- colnames(d_obs_resid) <- ids
          d_obs_resid_long <- d_obs_resid %>%
            reshape2::melt(value.name = paste0(var_name, "_dist_resid")) %>%
            dplyr::rename(id1 = "Var1", id2 = "Var2")
          # combine
          dplyr::full_join(
            d_obs_long, d_obs_resid_long, by = c("id1", "id2")
          ) %>% tibble::as_tibble()
        } else {
          d_obs_long %>% tibble::as_tibble()
        }
      }
    )
  purrr::reduce(
    d_obs_long_list,
    function(x, y) {
      dplyr::full_join(
        x, y, by = c("id1", "id2")
      )
    }
  )
}

#' @rdname pairwise_distances
#' @export
bin_pairwise_distances <- function(
    pairwise_distances,
    geo_bin = 100,
    time_bin = 100,
    per_bin_operation = function(x) { 0.5*mean(x^2, na.rm = T) }
) {
  # input check
  checkmate::assert_class(pairwise_distances, "mobest_pairwisedistances")
  # perform binning
  pairwise_distances %>%
    dplyr::mutate(
      geo_dist_cut = (cut(
        .data[["geo_dist"]],
        breaks = unique(c(seq(0, max(.data[["geo_dist"]]), geo_bin), max(.data[["geo_dist"]]))),
        include.lowest	= T, labels = F
      ) * geo_bin) - geo_bin/2,
      time_dist_cut = (cut(
        .data[["time_dist"]],
        breaks = unique(c(seq(0, max(.data[["time_dist"]]), time_bin), max(.data[["time_dist"]]))),
        include.lowest	= T, labels = F
      ) * time_bin) - time_bin/2
    ) %>%
    dplyr::group_by(.data[["geo_dist_cut"]], .data[["time_dist_cut"]]) %>%
    dplyr::summarise(
      dplyr::across(
        -tidyselect::any_of(c("id1", "id2", "geo_dist", "time_dist")),
        per_bin_operation
      ),
      n = dplyr::n(),
      .groups	= "drop"
    ) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_empiricalvariogram")
}
