#' Input data type constructors
#'
#' Functions to create the main input data types for the mobest package.
#' See the README for example code.
#'
#' @param id Vector. IDs of the observation points
#' @param x Numeric vector. Spatial x-axis coordinates
#' @param y Numeric vector. Spatial y-axis coordinates
#' @param z Numeric vector. Temporal positions
#' @param dsx Double. Kernel lengthscale parameter for the x dimension (spatial x-axis). See \code{?laGP::newGP} for more info
#' @param dsy Double. Kernel lengthscale parameter for the y dimension (spatial y-axis)
#' @param dt Double. Kernel lengthscale parameter for the z dimension (temporal axis)
#' @param g Double. Kernel nugget parameter
#' @param on_residuals Logical. In the field calculation down the pipeline: Should a linear model take out the main trends before the kriging interpolation?
#' @param auto Logical. In the field calculation down the pipeline:
#' Should the lengthscale and nugget values be automatically determined by laGPs
#' maximum likelihood algorithm? See \code{?laGP::mleGPsep} for more info
#' @param ... Different inputs (see examples in README)
#' @param .names Vector. Names of different object iterations
#'
#' @return Different data types for specific applications.
#'
#' @rdname input_data_constructors
#' @export
create_obs <- function(..., .names = NULL) {
  obs <- list(...)
  if (!is.null(.names)) { names(obs) <- .names }
  # check list
  checkmate::assert_list(obs, types = "numeric", names = "strict")
  checkmate::assert_true(
    purrr::map_int(obs, length) %>% unique %>% length %>% magrittr::equals(1)
  )
  # compile tibble
  dplyr::bind_cols(obs) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_observations")
}

#' @rdname input_data_constructors
#' @export
create_geopos <- function(id, x, y, ...) {
  # input check
  checkmate::assert_atomic_vector(id, any.missing = F, unique = T)
  checkmate::assert_numeric(x)
  checkmate::assert_numeric(y)
  checkmate::assert_true(
    purrr::map_int(list(id, x, y, ...), length) %>% unique %>% length %>% magrittr::equals(1)
  )
  # compile tibble
  tibble::tibble(id = id, x = x, y = y, ...) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_spatialpositions")
}

#' @rdname input_data_constructors
#' @export
create_geopos_multi <- function(..., .names = NULL) {
  geopos <- list(...)
  if (!is.null(.names)) { names(geopos) <- .names }
  # input check
  checkmate::assert_list(geopos, types = "mobest_spatialpositions", names = "strict")
  checkmate::assert_true(
    purrr::map_int(geopos, nrow) %>% unique %>% length %>% magrittr::equals(1)
  )
  checkmate::assert_true(
    purrr::map_lgl(geopos, function(x) { all(x[["id"]] == geopos[[1]]$id) }) %>% all()
  )
  # compile output data structure
  geopos
}

geopos_to_spatpos <- function(geopos, z) {
  # input check
  checkmate::assert_class(geopos, classes = "mobest_spatialpositions")
  checkmate::assert_numeric(z)
  # expand grid
  spatpos <- geopos %>%
    tidyr::crossing(tibble::tibble(z = z))
  # compile output
  mobest::create_spatpos(
    id = 1:nrow(spatpos),
    x = spatpos$x,
    y = spatpos$y,
    z = spatpos$z,
    geo_id = spatpos$id
  ) %>%
    dplyr::bind_cols(spatpos %>% dplyr::select(-id, -x, -y, -z))
}

#' @rdname input_data_constructors
#' @export
create_spatpos <- function(id, x, y, z, ...) {
  # input check
  checkmate::assert_atomic_vector(id, any.missing = F, unique = T)
  checkmate::assert_numeric(x)
  checkmate::assert_numeric(y)
  checkmate::assert_numeric(z)
  checkmate::assert_true(
    purrr::map_int(list(id, x, y, z, ...), length) %>% unique %>% length %>% magrittr::equals(1)
  )
  # compile tibble
  tibble::tibble(id = id, x = x, y = y, z = z, ...) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_spatiotemporalpositions")
}

#' @rdname input_data_constructors
#' @export
create_spatpos_multi <- function(..., .names = NULL) {
  spatpos <- list(...)
  if (!is.null(.names)) { names(spatpos) <- .names }
  # input check
  checkmate::assert_list(spatpos, types = "mobest_spatiotemporalpositions", names = "strict")
  checkmate::assert_true(
    purrr::map_int(spatpos, nrow) %>% unique %>% length %>% magrittr::equals(1)
  )
  checkmate::assert_true(
    purrr::map_lgl(spatpos, function(x) { all(x[["id"]] == spatpos[[1]]$id) }) %>% all()
  )
  # compile output data structure
  spatpos
}

#' @rdname input_data_constructors
#' @export
create_kernel <- function(dsx, dsy, dt, g, on_residuals = T, auto = F) {
  # input check
  checkmate::assert_number(dsx, lower = 0)
  checkmate::assert_number(dsy, lower = 0)
  checkmate::assert_number(dt, lower = 0)
  checkmate::assert_number(g, lower = 0)
  checkmate::assert_logical(on_residuals, len = 1)
  checkmate::assert_logical(auto, len = 1)
  # compile output data structure
  list(
    dsx = dsx,
    dsy = dsy,
    dt = dt,
    g = g,
    on_residuals = on_residuals,
    auto = auto
  ) %>%
    magrittr::set_class("mobest_kernel")
}

#' @rdname input_data_constructors
#' @export
create_kernset <- function(..., .names = NULL) {
  kernels <- list(...)
  if (!is.null(.names)) { names(kernels) <- .names }
  # input check
  checkmate::assert_list(kernels, types = "mobest_kernel", names = "strict")
  # compile output data structure
  kernels %>% magrittr::set_class("mobest_kernelsetting")
}

#' @rdname input_data_constructors
#' @export
create_kernset_multi <- function(..., .names = NULL) {
  kernels_multi <- list(...)
  if (!is.null(.names)) { names(kernels_multi) <- .names }
  # input check
  checkmate::assert_list(kernels_multi, types = "mobest_kernelsetting", names = "strict")
  # compile output data structure
  kernels_multi
}
