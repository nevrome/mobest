#' Input data type constructors
#'
#' Functions to create the main input data types for the mobest package.
#'
#' @param id Vector. Identifiers of the individual data points
#' @param x Numeric vector. Spatial x-axis coordinates
#' @param y Numeric vector. Spatial y-axis coordinates
#' @param z Numeric vector. Temporal positions
#' @param area An object of class \code{sf}. Polygons where the spatial grid should
#' be constructed
#' @param spatial_cell_size Numeric. Size of the output spatial grid cells in the unit of
#' \code{area}. See \code{?sf::st_make_grid} for more info
#' @param dsx Double. Kernel lengthscale parameter for the x dimension (spatial x-axis). See \code{?laGP::newGP} for more information
#' @param dsy Double. Kernel lengthscale parameter for the y dimension (spatial y-axis)
#' @param dt Double. Kernel lengthscale parameter for the z dimension (temporal axis)
#' @param g Double. Kernel nugget parameter
#' @param on_residuals Logical. In the field calculation:
#' Should a linear model take over the main trends before the kriging interpolation?
#' @param auto Logical. In the field calculation:
#' Should the lengthscale and nugget values be automatically determined by laGPs
#' maximum likelihood algorithm? See \code{?laGP::mleGPsep} for more info
#' @param ... Different inputs - usually individual elements to be merged in a list type
#' @param .names Vector. Names of different object iterations
#'
#' @name data_types
NULL

#### spatial coordinates ####

#' @rdname data_types
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

#' @rdname data_types
#' @export
create_prediction_grid <- function(area, spatial_cell_size) {
  # dependency check
  check_if_packages_are_available("sf")
  # input checks
  checkmate::assert_class(area, classes = "sf")
  # prepare grid
  sf::st_agr(area) <- "constant"
  space_grid_rect_sf <- area %>%
    sf::st_make_grid(cellsize = spatial_cell_size, what = "centers") %>%
    sf::st_sf()
  sf::st_agr(space_grid_rect_sf) <- "constant"
  space_grid_sf <- sf::st_intersection(space_grid_rect_sf, area)
  space_grid <- space_grid_sf %>%
    dplyr::mutate(
      x = sf::st_coordinates(.)[,1],
      y = sf::st_coordinates(.)[,2]
    ) %>%
    sf::st_drop_geometry() %>%
    dplyr::select(.data[["x"]], .data[["y"]])
  # compile output
  mobest::create_geopos(
    id = 1:nrow(space_grid),
    x = space_grid$x,
    y = space_grid$y
  )
}

#' @rdname data_types
#' @export
create_geopos_multi <- function(..., .names = NULL) {
  tibble_multi_function_factory(
    "mobest_spatialpositions",
    "mobest_spatialpositions_multi",
    T,T
  )(..., .names = .names)
}

#### spatiotemporal coordinates ####


#' @param geopos An object of class `mobest_spatialpositions`.
#'
#' @rdname data_types
#' @export
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
    dplyr::bind_cols(
      spatpos %>% dplyr::select(
        -.data[["id"]], -.data[["x"]], -.data[["y"]], -.data[["z"]]
      )
    )
}

#' @rdname data_types
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

#' @rdname data_types
#' @export
create_spatpos_multi <- function(..., .names = NULL) {
  tibble_multi_function_factory(
    "mobest_spatiotemporalpositions",
    "mobest_spatiotemporalpositions_multi",
    T,T
  )(..., .names = .names)
}

#### genetic coordinates ####

#' @rdname data_types
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

#' @rdname data_types
#' @export
create_obserror <- function(..., .names = NULL) {
  obs <- list(...)
  if (!is.null(.names)) { names(obs) <- .names }
  # check list
  checkmate::assert_list(obs, types = "numeric", names = "strict")
  checkmate::assert_true(all(grepl("_sd", names(obs))))
  checkmate::assert_true(
    purrr::map_int(obs, length) %>% unique %>% length %>% magrittr::equals(1)
  )
  # compile tibble
  dplyr::bind_cols(obs) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_observations_error")
}

#' @param obs An object of class `mobest_observations`.
#' @param obserror An object of class `mobest_observations_error`.
#'
#' @rdname data_types
#' @export
create_obs_obserror <- function(obs, obserror) {
  # input check
  checkmate::assert_class(obs, "mobest_observations")
  checkmate::assert_class(obserror, "mobest_observations_error")
  checkmate::assert_true(
    all(names(obserror) == paste0(names(obs), "_sd"))
  )
  checkmate::assert_true(nrow(obs) == nrow(obserror))
  # compile tibble
  dplyr::bind_cols(obs, obserror) %>%
    tibble::new_tibble(., nrow = nrow(.), class = "mobest_observationswitherror")
}

#' @rdname data_types
#' @export
create_obs_multi <- function(..., .names = NULL) {
  tibble_multi_function_factory(
    "mobest_observations",
    "mobest_observations_multi",
    T,F
  )(..., .names = .names)
}

#' @rdname data_types
#' @export
create_obs_obserror_multi <- function(..., .names = NULL) {
  tibble_multi_function_factory(
    "mobest_observationswitherror",
    "mobest_observationswitherror_multi",
    T,F
  )(..., .names = .names)
}

#### kernel settings ####

#' @rdname data_types
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

#' @rdname data_types
#' @export
create_kernset <- function(..., .names = NULL) {
  tibble_multi_function_factory(
    "mobest_kernel",
    "mobest_kernelsetting",
    F,F
  )(..., .names = .names)
}

#' @rdname data_types
#' @export
create_kernset_multi <- function(..., .names = NULL) {
  tibble_multi_function_factory(
    "mobest_kernelsetting",
    "mobest_kernelsetting_multi",
    F,F
  )(..., .names = .names)
}

#### helper functions ####

# this function produces other constructor functions
tibble_multi_function_factory <- function(single_type, multi_type, is_df = F, has_id = F) {
  function(..., .names = NULL) {
    multi <- list(...)
    if (!is.null(.names)) { names(multi) <- .names }
    # input check
    checkmate::assert_list(
      multi, types = single_type,
      any.missing = F, min.len = 1,
      names = "strict"
    )
    if (is_df) {
      checkmate::assert_true(
        purrr::map_int(multi, nrow) %>% unique %>% length %>% magrittr::equals(1)
      )
    }
    if (has_id) {
      checkmate::assert_true(
        purrr::map_lgl(multi, function(x) { all(x[["id"]] == multi[[1]]$id) }) %>% all()
      )
    }
    # compile output data structure
    multi %>% magrittr::set_class(c(multi_type, "list"))
  }
}
