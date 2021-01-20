#' Input data type constructors
#'
#' Functions to create the main input data types for the mobest package.
#' See the README for example code.
#'
#' @param id vector. IDs of the observation points
#' @param x numeric vector. Spatial x-axis coordinates
#' @param y numeric vector. Spatial y-axis coordinates
#' @param z numeric vector. Temporal positions
#' @param d 3 element numeric vector. Kernel lengthscale parameter for the x, y and z dimension
#' @param ds numeric vector. Different kernel lengthscale parameters for the x and y dimension
#' (isotropic)
#' @param dt numeric vector. Different kernel lengthscale parameters for the z dimension
#' @param g numeric vector. Kernel nugget parameter
#' @param on_residuals logical. Should the field calculated with these settings be calculated
#' on real values or the residuals of a linear model?
#' @param auto logical. Should the field calculated with these settings be calculated
#' on exactly these values for d and g, or should laGP attempt to estimate them from the data.
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
create_kernset <- function(d, g, on_residuals = T, auto = F, ...) {
  # input check
  checkmate::assert_numeric(d, len = 3)
  checkmate::assert_number(g)
  checkmate::assert_logical(on_residuals, len = 1)
  checkmate::assert_logical(auto, len = 1)
  # compile output data structure
  list(
    d = d,
    g = g,
    on_residuals = on_residuals,
    auto = auto,
    ...
  ) %>%
    magrittr::set_class("mobest_kernelsetting")
}

create_kernset_multi <- function(d, g, on_residuals = T, auto = F, ..., it) {
  # input check
  checkmate::assert_list(d)
  checkmate::assert_vector(g)
  checkmate::assert_logical(on_residuals)
  checkmate::assert_logical(auto)
  checkmate::assert_atomic_vector(it, any.missing = F, unique = T)
  if (list(d, g, on_residuals, auto, ...) %>%
      purrr::some(function(x) { length(x) != length(it) && length(x) != 1 }))
  { stop("Each input list must have identical length") }
  # compile output data structure
  purrr::pmap(list(d, g, on_residuals, auto, ...), mobest::create_kernset) %>%
    magrittr::set_names(it)
}

#' @rdname input_data_constructors
#' @export
create_kernset_cross <- function(ds, dt, g, on_residuals = T, auto = F) {
  # input check
  checkmate::assert_numeric(ds)
  checkmate::assert_numeric(dt)
  checkmate::assert_numeric(g)
  checkmate::assert_logical(on_residuals)
  checkmate::assert_logical(auto)
  # cross
  ks <- expand.grid(ds = ds, dt = dt, g = g)
  # create kernset list
  create_kernset_multi(
    d = purrr::map2(ks$ds, ks$dt, function(x, y) { c(x, x, y) }),
    g = as.list(ks$g),
    on_residuals = on_residuals,
    auto = auto,
    it = purrr::pmap_chr(
      list(ks$ds, ks$dt, ks$g),
      function(x, y, z) { paste("kernel", ps(x), ps(y), ps(z), sep = "_") }
    )
  )
}

ps <- function(x) {
  suppressWarnings(format(x, scientific = FALSE, decimal.mark = ""))
}

