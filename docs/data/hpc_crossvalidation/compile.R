kernel_grid <- purrr::map_dfr(
  list.files(
    path = "docs/data/hpc_crossvalidation",
    pattern = "*.csv",
    full.names = TRUE
  ),
  function(x) {
    readr::read_csv(x, show_col_types = FALSE)
  }
)
