### Temporal resampling as a permutation application

#### Drawing ages for each sample

library(magrittr)
library(ggplot2)

samples_advanced <- readr::read_csv("docs/data/samples_advanced.csv")

radiocarbon_date_sumcal <- function(ages, sds, cal_curve) {
  bol <- 1950
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
          radiocarbon_date_sumcal(
            ages = as.numeric(strsplit(c14bps, ";")[[1]]),
            sds = as.numeric(strsplit(c14sds, ";")[[1]]),
            cal_curve = "intcal20"
          )
        } else {
          contextual_date_uniform(
            startbcad = context_start,
            stopbcad = context_stop
          )
        }
      }
    )
  )

age_resampling_runs <- 2

set.seed(123)

samples_with_age_samples <- samples_with_age_densities %>%
  dplyr::mutate(
    Date_BC_AD_Samples = purrr::map(
      Date_BC_AD_Prob, function(x) {
          sample(x = x$age, size = age_resampling_runs, prob = x$sum_dens, replace = T)
      }
    )
  )

#### Applying the similarity search

load("docs/data/simple_objects_snapshot.RData")

dep_multi <- mobest::create_obs_multi(d = dep)
kernset_multi <- mobest::create_kernset_multi(k = kernset)
search_dep_multi <- mobest::create_obs_multi(d = search_dep)

ind_multi <- do.call(
  mobest::create_spatpos_multi,
  c(
    purrr::map(
      seq_len(age_resampling_runs), function(age_resampling_run) {
        mobest::create_spatpos(
          id = samples_with_age_samples$Sample_ID,
          x  = samples_with_age_samples$x,
          y  = samples_with_age_samples$y,
          z  = purrr::map_int(
            samples_with_age_samples$Date_BC_AD_Samples,
            function(x) {x[age_resampling_run]}
          )
        )
      }
    ),
    list(.names = paste0("age_resampling_run_", seq_len(age_resampling_runs)))
  )
)

search_samples <- samples_with_age_samples %>%
  dplyr::filter(
    Sample_ID == "Stuttgart_published.DG"
  )

search_samples_multi <- do.call(
  mobest::create_spatpos_multi,
  c(
    purrr::map(
      seq_len(age_resampling_runs), function(age_resampling_run) {
        mobest::create_spatpos(
          id = search_samples$Sample_ID,
          x  = search_samples$x,
          y  = search_samples$y,
          z  = search_samples$Date_BC_AD_Median
        )
      }
    ),
    list(.names = paste0("age_resampling_run_", seq_len(age_resampling_runs)))
  )
)

search_result <- mobest::locate_multi(
  independent        = ind_multi,
  dependent          = dep_multi,
  kernel             = kernset_multi,
  search_independent = search_samples_multi,
  search_dependent   = search_dep_multi,
  search_space_grid  = spatial_pred_grid,
  search_time        = -6800,
  search_time_mode   = "absolute"
)

search_product <- mobest::multiply_dependent_probabilities(search_result)

p <- ggplot() +
  geom_raster(
    data = search_product,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  scale_fill_viridis_c() +
  geom_sf(
    data = research_area_3035,
    fill = NA, colour = "red",
    linetype = "solid", linewidth = 1
  ) +
  geom_point(
    data = search_samples,
    mapping = aes(x, y),
    colour = "red"
  ) +
  ggtitle(
    label = "<Stuttgart> ~5250BC",
    subtitle = "Early Neolithic (Linear Pottery Culture) - Lazaridis et al. 2014"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  ) +
  guides(
    fill = guide_colourbar(title = "Similarity\nsearch\nprobability")
  ) +
  annotate(
    "text",
    label = "6800BC",
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.5
  ) +
  facet_wrap(~independent_table_id)

# ggsave(
#   filename = "docs/img/temporal_resampling/search_map_two_resampling_runs.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

search_sum <- mobest::fold_probabilities_per_group(search_product)

p <- ggplot() +
  geom_raster(
    data = search_sum,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  scale_fill_viridis_c() +
  geom_sf(
    data = research_area_3035,
    fill = NA, colour = "red",
    linetype = "solid", linewidth = 1
  ) +
  geom_point(
    data = search_samples,
    mapping = aes(x, y),
    colour = "red"
  ) +
  ggtitle(
    label = "<Stuttgart> ~5250BC",
    subtitle = "Early Neolithic (Linear Pottery Culture)\nLazaridis et al. 2014"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  ) +
  guides(
    fill = guide_colourbar(title = "Similarity\nsearch\nprobability")
  ) +
  annotate(
    "text",
    label = "6800BC",
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.5
  )

# ggsave(
#   filename = "docs/img/temporal_resampling/search_map_combined_resampling_runs.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )
