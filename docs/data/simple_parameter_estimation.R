library(ggplot2)
library(magrittr)

# Parameter estimation for optimal ancestry interpolation

## Using a subset of the variogram to estimate the nugget parameter

samples_projected <- readr::read_csv("docs/data/samples_projected.csv")

### Determining pairwise distances

distances_all <- mobest::calculate_pairwise_distances(
  independent = mobest::create_spatpos(
    id = samples_projected$Sample_ID,
    x = samples_projected$x,
    y = samples_projected$y,
    z = samples_projected$Date_BC_AD_Median
  ),
  dependent = mobest::create_obs(
    C1 = samples_projected$MDS_C1,
    C2 = samples_projected$MDS_C2
  )
)

p1 <- ggplot() +
  geom_bin2d(
    data = distances_all,
    mapping = aes(x = geo_dist, y = obs_dist_total),
    bins = 30
  ) +
  scale_fill_viridis_c() +
  theme_bw()

p2 <- ggplot() +
  geom_bin2d(
    data = distances_all,
    mapping = aes(x = time_dist, y = obs_dist_total),
    bins = 30
  ) +
  scale_fill_viridis_c() +
  theme_bw()

p <- cowplot::plot_grid(p1, p2)

# ggsave(
#   filename = "docs/img/estimation/distance_correlation.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

### Summarizing distances in an empirical variogram

variogram <- mobest::bin_pairwise_distances(
  distances_all,
  geo_bin = 100, time_bin = 100
)

### Estimating the nugget parameter

distances_for_nugget <- distances_all %>%
  # remove auto-distances
  dplyr::filter(id1 != id2) %>%
  # filter for small temporal and spatial pairwise distances
  dplyr::filter(time_dist < 50 & geo_dist < 50) %>%
  # transform the residual dependent variable distances
  # into a long format table
  tidyr::pivot_longer(
    cols = tidyselect::ends_with("_resid"),
    names_to = "dist_type", values_to = "dist_val"
  ) %>%
  # rescale the distances to relative proportions
  dplyr::mutate(
    dist_val_adjusted = dplyr::case_when(
      dist_type == "C1_dist_resid" ~
        0.5*(dist_val^2 / stats::var(samples_projected$MDS_C1)),
      dist_type == "C2_dist_resid" ~
        0.5*(dist_val^2 / stats::var(samples_projected$MDS_C2))
    )
  )

estimated_nuggets <- distances_for_nugget %>%
  dplyr::group_by(dist_type) %>%
  dplyr::summarise(nugget = mean(dist_val_adjusted, na.rm = T))

p <- ggplot() +
  geom_violin(
    data = distances_for_nugget,
    mapping = aes(x = dist_type, y = dist_val_adjusted, fill = dist_type),
    linewidth = 0.5,
    width = 0.8
  ) +
  geom_boxplot(
    data = distances_for_nugget,
    mapping = aes(x = dist_type, y = dist_val_adjusted),
    width = 0.1, outlier.size = 1
  ) +
  geom_point(
    data = estimated_nuggets,
    mapping = aes(x = dist_type, y = nugget),
    size = 4, shape = 18
  ) +
  geom_point(
    data = estimated_nuggets,
    mapping = aes(x = dist_type, y = nugget),
    size = 6, shape = "|"
  ) +
  geom_text(
    data = estimated_nuggets,
    mapping = aes(x = dist_type, y = nugget, label = paste0("mean: ~", round(nugget, 3))),
    nudge_x = -0.5
  ) +
  coord_flip() +
  theme_bw() +
  guides(fill = "none") +
  xlab("ancestry component distance type") +
  ylab("pairwise half mean squared normalized residual distance") +
  scale_y_log10(labels = scales::comma) +
  scale_x_discrete(limits = rev(unique(distances_for_nugget$dist_type)))

# ggsave(
#   filename = "docs/img/estimation/nuggets.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

## Finding optimal lengthscale parameters with crossvalidation

### Basic setup

### Analyzing the crossvalidation results

### HPC setup for large lengthscale parameter spaces
