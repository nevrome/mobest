#' Title
#'
#' @param ...
#'
#' @return
#'
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

merge_spatpos_obs <- function(spatpos, obs) {
  dplyr::inner_join(
    spatpos, obs, by = "id"
  )
}

mean_spatpos_uncertain <- function(spatpos_uncertain) {
  # input check
  checkmate::assert_class(spatpos_uncertain, "mobest_uncertain_spatiotemporalpositions")
  checkmate::assert_true(length(spatpos_uncertain) > 1)
  # calculate mean
  purrr::reduce(
    spatpos_uncertain, function(a, b) {
      create_spatpos(
        id = a$id,
        x = a$x+b$x,
        y = a$y+b$y,
        z = a$z+b$z
      )
    }
  ) %>% dplyr::mutate(
    dplyr::across(c("x", "y", "z"), function(x) { x/length(spatpos_uncertain) })
  )
}

#' Title
#'
#' @param id
#' @param x
#' @param y
#' @param z
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param id
#' @param x
#' @param y
#' @param z
#' @param name
#'
#' @return
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

#' Title
#'
#' @param d
#' @param g
#' @param on_residuals
#' @param auto
#'
#' @return
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

#' Title
#'
#' @param d
#' @param g
#' @param ...
#' @param it
#'
#' @return
#' @export
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

#' create_kernel_grid
#'
#' @param ds test
#' @param dt test
#' @param g test
#'
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

