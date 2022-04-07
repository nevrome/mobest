#' Input data type constructors
#'
#' Functions to create the main input data types for the mobest package.
#' See the README for example code.
#'
#' @param id Vector. IDs of the observation points
#' @param x Numeric vector. Spatial x-axis coordinates
#' @param y Numeric vector. Spatial y-axis coordinates
#' @param z Numeric vector. Temporal positions
#' @param d 3 element numeric vector. Kernel lengthscale parameter for the x, y and z dimension. See \code{?laGP::newGP} for more info
#' @param ds numeric vector. Different kernel lengthscale parameters for the x and y dimension
#' (isotropic)
#' @param dt numeric vector. Different kernel lengthscale parameters for the z dimension
#' @param g numeric vector. Kernel nugget parameter
#' @param on_residuals logical. In the field calculation down the pipeline: Should a linear model take out the main trends before the kriging interpolation?
#' @param auto logical. In the field calculation down the pipeline:
#' Should the lengthscale and nugget values be automatically determined by laGPs
#' maximum likelihood algorithm? See \code{?laGP::mleGPsep} for more info
#' @param ... vector. Other settings to be added to the output object
#' @param it vector. Names of the different object iterations
#'
#' @return Different data types for specific applications.
#'
#' @rdname input_data_constructors
#' @export
create_obs <- function(...) {
  # prepare list
  res <- list(...)
  # check list
  checkmate::assert_names(names(res), type = "strict")
  purrr::walk(res, checkmate::assert_numeric)
  if (res %>%
      purrr::some(function(x) { length(x) != length(res[[1]]) }))
  { stop("Each input vector must have identical length") }
  # return list
  class(res) <- c("mobest_observations", class(res))
  return(res)
}

#' @rdname input_data_constructors
#' @export
create_spatpos <- function(id, x, y, z, ...) {
  # input check
  checkmate::assert_atomic_vector(id, any.missing = F, unique = T)
  checkmate::assert_numeric(x)
  checkmate::assert_numeric(y)
  checkmate::assert_numeric(z)
  if (list(id, x, y, z, ...) %>%
      purrr::some(function(x) { length(x) != length(id) && length(x) != 1 }))
  { stop("Each vector input vector must be of identical length") }
  # compile tibble
  tibble::tibble(id = id, x = x, y = y, z = z, ...) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_spatiotemporalpositions")
}

#' @rdname input_data_constructors
#' @export
create_spatpos_multi <- function(id, x, y, z, it) {
  # input check
  checkmate::assert_vector(id)
  checkmate::assert_list(x)
  checkmate::assert_list(y)
  checkmate::assert_list(z)
  checkmate::assert_atomic_vector(it, any.missing = F, unique = T)
  if (list(x, y, z) %>%
      purrr::some(function(x) { length(x) != length(it) }))
    { stop("Each input list must have identical length") }
  # compile list of tibbles
  list(
    id = rep(list(id), length(x)),
    x = x, y = y, z = z
  ) %>%
    purrr::pmap(
      function(id, x, y, z, it) {
        create_spatpos(id = id, x = x, y = y, z = z)
      }
    ) %>%
    magrittr::set_names(it)
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
create_kernset <- function(...) {
  kernels <- list(...)
  # input check
  checkmate::assert_list(kernels, types = "mobest_kernel", names = "strict")
  # compile output data structure
  kernels %>% magrittr::set_class("mobest_kernelsetting")
}

#' @rdname input_data_constructors
#' @export
create_kernset_multi <- function(...) {
  kernels_multi <- list(...)
  # input check
  checkmate::assert_list(kernels_multi, types = "mobest_kernelsetting", names = "strict")
  # compile output data structure
  kernels_multi
}

ps <- function(x) {
  suppressWarnings(format(x, scientific = FALSE, decimal.mark = ""))
}

