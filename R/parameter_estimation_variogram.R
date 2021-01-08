#' Title
#'
#' @param independent
#' @param dependent
#' @param m_to_km
#'
#' @return
#'
#' @export
calculate_pairwise_distances <- function(independent, dependent, m_to_km = T) {
  # input check and transformation
  checkmate::assert_class(independent, "mobest_spatiotemporalpositions")
  checkmate::assert_class(dependent, "mobest_observations")
  dependent <- dplyr::bind_cols(dependent)
  if (nrow(independent) != nrow(dependent)) {
    stop("independent and dependent must have the same number of rows")
  }
  ids <- independent$id
  # geo distance
  d_geo <- dist(independent %>% dplyr::select(.data[["x"]], .data[["y"]]), "euclidean") %>% as.matrix()
  rownames(d_geo) <- colnames(d_geo) <- ids
  d_geo_long <- d_geo %>%
    reshape2::melt(value.name = "geo_dist") %>%
    dplyr::mutate(
      # m to km
      geo_dist = if (m_to_km) {.data[["geo_dist"]]/1000} else {.data[["geo_dist"]]}
    )
  # time distance
  d_time <- dist(independent %>% dplyr::select(.data[["z"]]), "euclidean") %>% as.matrix()
  rownames(d_time) <- colnames(d_time) <- ids
  d_time_long <- d_time %>% reshape2::melt(value.name = "time_dist")
  # obs total distance
  d_obs_total <- dist(dependent, "euclidean") %>% as.matrix()
  rownames(d_obs_total) <- colnames(d_obs_total) <- ids
  d_obs_total_long <- d_obs_total %>% reshape2::melt(value.name = "obs_dist_total")
  # obs individual distance
  var_names <- names(dependent)
  d_obs_long_list <- var_names %>%
    purrr::map(function(var_name) {
      # d_obs
      d_obs <- dist(dependent[[var_name]], "euclidean") %>% as.matrix()
      rownames(d_obs) <- colnames(d_obs) <- ids
      d_obs_long <- d_obs %>% reshape2::melt(value.name = paste0(var_name, "_dist"))
      # d_obs_resid
      model <- lm(dependent[[var_name]] ~ independent$x + independent$y + independent$z)
      model_residuals <- residuals(model)
      d_obs_resid <- as.matrix(dist(model_residuals, "euclidean"))
      rownames(d_obs_resid) <- colnames(d_obs_resid) <- ids
      d_obs_resid_long <- d_obs_resid %>% reshape2::melt(value.name = paste0(var_name, "_dist_resid"))
      # combine
      dplyr::full_join(
        d_obs_long, d_obs_resid_long, by = c("Var1", "Var2")
      )
    })
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

#' Title
#'
#' @param x
#' @param geo_bin
#' @param time_bin
#'
#' @return
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
