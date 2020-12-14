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
create_spatpos <- function(id, x, y, z) {
  # input check
  checkmate::assert_atomic_vector(id, any.missing = F, unique = T)
  checkmate::assert_numeric(x)
  checkmate::assert_numeric(y)
  checkmate::assert_numeric(z)
  if (list(id, x, y, z) %>%
      purrr::some(function(x) { length(x) != length(id) }))
  { stop("Each vector in id, x, y and z must have identical length") }
  # compile tibble
  tibble::tibble(id = id, x = x, y = y, z = z) %>%
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
    { stop("x, y, z, and it must have identical length") }
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

get_var_names <- function(obs) {
  t_obs <- colnames(obs)
  t_obs[t_obs != "id"]
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
create_kernset <- function(d, g, on_residuals, auto) {
  list(
    d = d,
    g = g,
    on_residuals = on_residuals,
    auto = auto
  ) %>%
    magrittr::set_class("mobest_kernel_setting")
}
