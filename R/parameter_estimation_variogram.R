#' Prepare an empirical semivariogram
#'
#' @param independent An object of class mobest_spatiotemporalpositions
#' @param dependent An object of class mobest_observations
#' @param m_to_km Logical. Should distances be transformed from m to km (x/1000)
#' @param x An object of class mobest_pairwisedistances
#' @param geo_bin Numeric. Width of the spatial bins
#' @param time_bin Numeric. Width of the temporal bins
#' @param with_resid Logical. Calculate distances also on the residuals of a spatiotemporal linear model
#'
#' @rdname variogram
#' @export
calculate_pairwise_distances <- function(independent, dependent, m_to_km = T, with_resid = T) {
  # input check and transformation
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  checkmate::assert_class(dependent, "mobest_observations")
  if (nrow(independent) != nrow(dependent)) {
    stop("independent and dependent must have the same number of rows")
  }
  ids <- independent$id
  # geo distance
  d_geo_long <- calculate_geo_pairwise_distances(ids, independent, m_to_km = m_to_km)
  # time distance
  d_time_long <- calculate_time_pairwise_distances(ids, independent)
  # obs total distance
  d_obs_total <- stats::dist(dependent, "euclidean") %>% as.matrix()
  rownames(d_obs_total) <- colnames(d_obs_total) <- ids
  d_obs_total_long <- d_obs_total %>% reshape2::melt(value.name = "obs_dist_total")
  # obs individual distance
  d_obs_long_list <- calculate_dependent_pairwise_distances(ids, dependent, with_resid = with_resid, independent)
  # join different distances
  purrr::reduce(
    c(list(d_geo_long, d_time_long, d_obs_total_long), d_obs_long_list),
    function(x, y) {
      dplyr::full_join(
        x, y, by = c("Var1", "Var2")
      )
    }
  ) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_pairwisedistances")
}

#' @param ids Character vector. Identifier for the observations
#'
#' @rdname variogram
#' @export
calculate_geo_pairwise_distances <- function(ids, independent, m_to_km = T) {
  d_geo <- stats::dist(independent %>% dplyr::select(.data[["x"]], .data[["y"]]), "euclidean") %>% as.matrix()
  rownames(d_geo) <- colnames(d_geo) <- ids
  d_geo %>%
    reshape2::melt(value.name = "geo_dist") %>%
    dplyr::mutate(
      # m to km
      geo_dist = if (m_to_km) {.data[["geo_dist"]]/1000} else {.data[["geo_dist"]]}
    )
}

#' @rdname variogram
#' @export
calculate_time_pairwise_distances <- function(ids, independent) {
  d_time <- stats::dist(independent %>% dplyr::select(.data[["z"]]), "euclidean") %>% as.matrix()
  rownames(d_time) <- colnames(d_time) <- ids
  d_time %>% reshape2::melt(value.name = "time_dist")
}

#' @rdname variogram
#' @export
calculate_dependent_pairwise_distances <- function(ids, dependent, with_resid = F, independent = NULL) {
  if (with_resid & is.null(independent)) { stop("If with_resid, then independent can not be NULL") }
  var_names <- names(dependent)
  var_names %>%
    purrr::map(
      function(var_name) {
        # d_obs
        d_obs <- stats::dist(dependent[[var_name]], "euclidean") %>% as.matrix()
        rownames(d_obs) <- colnames(d_obs) <- ids
        d_obs_long <- d_obs %>% reshape2::melt(value.name = paste0(var_name, "_dist"))
        # d_obs_resid
        if (with_resid) {
          model <- stats::lm(dependent[[var_name]] ~ independent$x + independent$y + independent$z)
          model_residuals <- stats::residuals(model)
          d_obs_resid <- as.matrix(stats::dist(model_residuals, "euclidean"))
          rownames(d_obs_resid) <- colnames(d_obs_resid) <- ids
          d_obs_resid_long <- d_obs_resid %>% reshape2::melt(value.name = paste0(var_name, "_dist_resid"))
          # combine
          dplyr::full_join(
            d_obs_long, d_obs_resid_long, by = c("Var1", "Var2")
          )
        } else {
          d_obs_long
        }
      }
    )
}

#' @rdname variogram
#' @export
bin_pairwise_distances <- function(x, geo_bin = 100, time_bin = 100) {
  # input check
  checkmate::assert_class(x, "mobest_pairwisedistances")
  # perform binning
  x %>%
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
        -tidyselect::any_of(c("Var1", "Var2", "geo_dist", "time_dist")),
        function(x) { 0.5*mean(x^2, na.rm = T) }
      ),
      n = dplyr::n(),
      .groups	= "drop"
    ) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_empiricalvariogram")
}
