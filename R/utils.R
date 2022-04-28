check_if_packages_are_available <- function(x) {
  if (
    x %>%
    sapply(function(x) {requireNamespace(x, quietly = TRUE)}) %>%
    all %>% `!`
  ) {
    stop(
      paste0(
        "R packages ",
        paste(x, collapse = ", "),
        " needed for this function to work. Please install with ",
        "install.packages(c('", paste(x, collapse = "', '"), "'))"
      ),
      call. = FALSE
    )
  }
}
