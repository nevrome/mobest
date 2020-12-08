create_spatpos <- function(id, x, y, z, name) {
  # input check and modifications
  checkmate::assert_vector(id)
  checkmate::assert_list(x)
  checkmate::assert_list(y)
  checkmate::assert_list(z)
  checkmate::assert_character(name)
  if (list(x, y, z) %>%
      purrr::some(function(x) { length(x) != length(name) }))
    { stop("id, x, y, z, and names must have identical length") }
  if (purrr::flatten(list(x, y, z)) %>%
      purrr::some(function(x) { length(x) != length(id) }))
    { stop("Each vector in id, x, y, z must have identical length") }
  # compile list of tibbles
  list(
    id = rep(list(id), length(x)),
    x = x, y = y, z = z
  ) %>%
    purrr::pmap(
      function(id, x, y, z) {
        tibble::tibble(
          id = id,
          x = x,
          y = y,
          z = z
        )
      }
    ) %>%
    magrittr::set_names(name) %>%
    magrittr::set_class("mobest_spatiotemporalposition")
}

create_obs <- function(id, ...) {
  # input check
  checkmate::assert_vector(id)
  if (list(...) %>%
      purrr::some(function(x) { length(x) != length(id) }))
    { stop("Each vector in ... must have the same length as id") }
  # compile obs tibble
  tibble::tibble(
    id = id
  ) %>%
    cbind(...) %>%
    magrittr::set_class("mobest_observation")
}

