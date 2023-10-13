library(magrittr)
library(ggplot2)

samples_advanced <- readr::read_csv("docs/data/samples_advanced.csv")

sumcal <- function(ages, sds, cal_curve) {
  bol <- 1950 # c14 reference zero
  raw_calibration_output <- Bchron::BchronCalibrate(
    ages      = ages,
    ageSds    = sds,
    calCurves = rep(cal_curve, length(ages))
  )
  density_tables <- purrr::map(
    raw_calibration_output,
    function(y) {
      tibble::tibble(
        age = as.integer(-y$ageGrid + bol),
        densities = y$densities
      )
    }
  )
  sum_density_table <- dplyr::bind_rows(density_tables) %>%
    dplyr::group_by(age) %>%
    dplyr::summarise(
      sum_dens = sum(densities)/length(density_tables),
      .groups = "drop"
    ) %>%
    dplyr::right_join(
      ., data.frame(age = min(.$age):max(.$age)), by = "age"
    ) %>%
    dplyr::mutate(sum_dens = tidyr::replace_na(sum_dens, 0)) %>%
    dplyr::arrange(age)
  return(sum_density_table)
}

contextual_date_uniform <- function(startbcad, stopbcad) {
  tibble::tibble(
    age = startbcad:stopbcad,
    sum_dens = 1/(length(startbcad:stopbcad))
  )
}

samples_with_age_densities <- samples_advanced %>%
  dplyr::mutate(
    Date_BC_AD_Prob = purrr::pmap(
      list(Date_BC_AD_Start, Date_BC_AD_Stop, C14_ages, C14_sds),
      function(context_start, context_stop, c14bps, c14sds) {
        if (!is.na(c14bps)) {
          sumcal(
            as.numeric(strsplit(c14bps, ";")[[1]]),
            as.numeric(strsplit(c14sds, ";")[[1]]),
            "intcal20"
          )
        } else {
          contextual_date_uniform(context_start, context_stop)
        }
      }
    )
  )

samples_with_age_samples <- samples_with_age_densities %>%
  dplyr::mutate(
    Date_BC_AD_Samples = purrr::map(
      Date_BC_AD_Prob, function(x) {
          sample(x = x$age, size = 100, prob = x$sum_dens, replace = T)
      }
    )
  )
