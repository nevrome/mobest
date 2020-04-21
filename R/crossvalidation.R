#' create_kernel_grid
#'
#' @param ds test
#' @param dt test
#' @param g test
#'
#' @export
create_kernel_grid <- function(ds, dt, g) {

  ks <- expand.grid(ds = ds, dt = dt, g = g)

  kernel_settings <- lapply(
      1:nrow(ks), function(i) {
        list(d = c(ks[["ds"]][i], ks[["ds"]][i], ks[["dt"]][i]), g = ks[["g"]][i], on_residuals = T, auto = F)
      }
    ) %>% setNames(
      sapply(
        1:nrow(ks), function(i) {
          paste0(ks[["ds"]][i]/1000, "_", ks[["dt"]][i], "_", ks[["g"]][i])
        }
      )
    )

  return(kernel_settings)

}
