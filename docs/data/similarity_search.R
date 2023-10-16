# A basic similarity search workflow

## Preparing the computational environment

library(magrittr)
library(ggplot2)

## Preparing the input data

### Generating the the spatial prediction grid

#### Defining the research area

research_area_4326 <- sf::st_polygon(
  list(
    cbind(
      c(35.91,11.73,-11.74,-15.47,37.06,49.26,49.56,35.91), # longitudes
      c(25.61,28.94, 31.77, 62.73,65.67,44.56,28.55,25.61)  # latitudes
    )
  )
) %>% sf::st_sfc(crs = 4326)

# mapview::mapview(research_area_4326)

worldwide_land_outline_4326 <- rnaturalearth::ne_download(
  scale = 50, type = 'land', category = 'physical',
  returnclass = "sf"
)

research_land_outline_4326 <- sf::st_intersection(
  worldwide_land_outline_4326,
  research_area_4326
)

p <- ggplot() + geom_sf(data = research_land_outline_4326)

# ggsave(
#   filename = "docs/img/basic/research_area_land_outline_4326.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

#### Projecting the spatial data

research_land_outline_3035 <- research_land_outline_4326 %>% sf::st_transform(crs = 3035)

p <- ggplot() + geom_sf(data = research_land_outline_3035)

# ggsave(
#   filename = "docs/img/basic/research_area_land_outline_3035.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

#### Creating the prediction grid

spatial_pred_grid <- mobest::create_prediction_grid(
  research_land_outline_3035,
  spatial_cell_size = 50000
)

p <- ggplot() +
  geom_sf(data = research_land_outline_3035) +
  geom_point(
    data = spatial_pred_grid,
    mapping = aes(x, y),
    color = "red",
    size = 0.25
  )

# ggsave(
#   filename = "docs/img/basic/spatial_prediction_grid.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

### Reading the input samples

samples_basic <- readr::read_csv("docs/data/samples_basic.csv")

samples_projected <- samples_basic %>%
  sf::st_as_sf(
    coords = c("Longitude", "Latitude"),
    crs = 4326
  ) %>%
  sf::st_transform(crs = 3035) %>%
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],
    y = sf::st_coordinates(.)[,2]
  ) %>%
  sf::st_drop_geometry()

p <- ggplot() +
  geom_sf(data = research_land_outline_3035) +
  geom_point(
    data = samples_projected,
    mapping = aes(x, y),
    color = "darkgreen",
    size = 0.25
  )

# ggsave(
#   filename = "docs/img/basic/samples_map.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

## Specifying the search sample

search_samples <- samples_projected %>%
  dplyr::filter(
    Sample_ID == "Stuttgart_published.DG"
  )

## Running mobest's interpolation and search function

### Building the input data for the interpolation

#### Independent and dependent positions

ind <- mobest::create_spatpos(
  id = samples_projected$Sample_ID,
  x  = samples_projected$x,
  y  = samples_projected$y,
  z  = samples_projected$Date_BC_AD_Median
)
dep <- mobest::create_obs(
  C1 = samples_projected$MDS_C1,
  C2 = samples_projected$MDS_C2
)

search_ind <- mobest::create_spatpos(
  id = search_samples$Sample_ID,
  x  = search_samples$x,
  y  = search_samples$y,
  z  = search_samples$Date_BC_AD_Median
)
search_dep <- mobest::create_obs(
  C1 = search_samples$MDS_C1,
  C2 = search_samples$MDS_C2
)

#### Kernel parameter settings

kernset <- mobest::create_kernset(
  C1 = mobest::create_kernel(
    dsx = 800 * 1000, dsy = 800 * 1000, dt = 800,
    g = 0.1
  ),
  C2 = mobest::create_kernel(
    dsx = 800 * 1000, dsy = 800 * 1000, dt = 800,
    g = 0.1
  )
)

#### Search positions

#### Normalization

### Calling mobest::locate()

search_result <- mobest::locate(
  independent        = ind,
  dependent          = dep,
  kernel             = kernset,
  search_independent = search_ind,
  search_dependent   = search_dep,
  search_space_grid  = spatial_pred_grid,
  search_time        = -6800,
  search_time_mode   = "absolute"
)

## Inspecting the computed results

### The mobest_locateoverview table

### Creating similarity probability maps for individual dependent variables

result_C1 <- search_result %>% dplyr::filter(dependent_var_id == "C1")

p_C1 <- ggplot() +
  geom_raster(
    data = result_C1,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  coord_fixed()

result_C2 <- search_result %>% dplyr::filter(dependent_var_id == "C2")

p_C2 <- ggplot() +
  geom_raster(
    data = result_C2,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  coord_fixed()

p <- cowplot::plot_grid(p_C1, p_C2, labels = c("C1", "C2"))

# ggsave(
#   filename = "docs/img/basic/search_map_simple_C1_C2.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

### Combining the information from multiple dependent variables

search_product <- mobest::multiply_dependent_probabilities(search_result)

p <- ggplot() +
  geom_raster(
    data = search_product,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  coord_fixed()

# ggsave(
#   filename = "docs/img/basic/search_map_simple_combined.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )



# Improving the similarity search map plot

spatial_pred_grid <- mobest::create_prediction_grid(
  research_land_outline_3035,
  spatial_cell_size = 20000
)

search_result <- mobest::locate(
  independent        = ind,
  dependent          = dep,
  kernel             = kernset,
  search_independent = search_ind,
  search_dependent   = search_dep,
  search_space_grid  = spatial_pred_grid,
  search_time        = -6800,
  search_time_mode   = "absolute"
)

search_product <- mobest::multiply_dependent_probabilities(search_result)

research_area_3035 <- research_area_4326 %>% sf::st_transform(3035)

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
    subtitle = "Early Neolithic (Linear Pottery Culture)\nLazaridis et al. 2014"
  ) +
  annotate(
    "text",
    label = "6800BC",
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.5
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  ) +
  guides(
    fill = guide_colourbar(title = "Similarity\nsearch\nprobability")
  )

# ggsave(
#   filename = "docs/img/basic/search_map_neat_combined.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

# save(
#   dep, kernset, spatial_pred_grid, search_dep, research_area_3035,
#   file = "docs/data/simple_objects_snapshot.RData"
# )



# Simple permutations

## Multiple search time slices

search_result <- mobest::locate(
  independent        = ind,
  dependent          = dep,
  kernel             = kernset,
  search_independent = search_ind,
  search_dependent   = search_dep,
  search_space_grid  = spatial_pred_grid,
  search_time        = c(-6800, -5700),
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
  facet_wrap(
    ~search_time,
    labeller = \(variable, value) {
      paste0("Search time: ", abs(value), "BC")
    }
  )

# ggsave(
#   filename = "docs/img/basic/search_map_two_timeslices.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

### Summarizing multiple search samples

origin_vectors <- mobest::determine_origin_vectors(search_product, search_time)

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
  geom_point(
    data = origin_vectors,
    mapping = aes(field_x, field_y),
    fill = "orange", shape = 24
  ) +
  geom_segment(
    data = origin_vectors,
    mapping = aes(
      x = search_x, y = search_y,
      xend = field_x, yend = field_y
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red"
  ) +
  geom_label(
    data = origin_vectors,
    mapping = aes(
      x = (field_x + search_x)/2, y = (field_y + search_y)/2,
      label = paste0(round(ov_dist/1000, -2), "km")
    ),
    fill = "white", colour = "red",
    size = 1.9, label.padding = unit(0.1, "lines")
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
  facet_wrap(
    ~search_time,
    labeller = \(variable, value) {
      paste0("Search time: ", abs(value), "BC")
    }
  )

# ggsave(
#   filename = "docs/img/basic/search_map_two_timeslices_with_ovs.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

## Multiple search samples

search_samples <- samples_projected %>%
  dplyr::filter(
    Sample_ID %in% c("Stuttgart_published.DG", "I5411")
  )

search_ind <- mobest::create_spatpos(
  id = search_samples$Sample_ID,
  x  = search_samples$x,
  y  = search_samples$y,
  z  = search_samples$Date_BC_AD_Median
)
search_dep <- mobest::create_obs(
  C1 = search_samples$MDS_C1,
  C2 = search_samples$MDS_C2
)

search_result <- mobest::locate(
  independent        = ind,
  dependent          = dep,
  kernel             = kernset,
  search_independent = search_ind,
  search_dependent   = search_dep,
  search_space_grid  = spatial_pred_grid,
  search_time        = -700,
  search_time_mode   = "relative"
)
search_product <- mobest::multiply_dependent_probabilities(search_result)
origin_vectors <- mobest::determine_origin_vectors(search_product)

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
    data = search_samples %>% dplyr::rename(search_id = Sample_ID),
    mapping = aes(x, y),
    colour = "red"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  ) +
  guides(
    fill = guide_colourbar(title = "Similarity\nsearch\nprobability")
  ) +
  facet_wrap(
    ~search_id,
    ncol = 2,
    labeller = labeller(
      search_id = c(
        "Stuttgart_published.DG" = paste(
          "<Stuttgart> ~5250BC",
          "Early Neolithic (Linear Pottery culture) - Lazaridis et al. 2014",
          "Search time: ~5950BC",
          sep = "\n"
        ),
        "I5411" = paste(
          "<I5411> ~6650BC",
          "Mesolithic (Iron Gates) - Mathieson et al. 2018",
          "Search time: ~7350BC",
          sep = "\n"
        )
      )
    )
  )

# ggsave(
#   filename = "docs/img/basic/search_map_two_samples.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )

p <- ggplot() +
  geom_sf(
    data = research_land_outline_3035,
    fill = "grey", color = NA
  ) +
  geom_sf(
    data = research_area_3035,
    fill = NA, colour = "red",
    linetype = "solid", linewidth = 1
  ) +
  geom_point(
    data = origin_vectors,
    mapping = aes(search_x, search_y),
    colour = "red"
  ) +
  geom_point(
    data = origin_vectors,
    mapping = aes(field_x, field_y),
    fill = "orange", shape = 24
  ) +
  geom_segment(
    data = origin_vectors,
    mapping = aes(
      x = search_x, y = search_y,
      xend = field_x, yend = field_y
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  )

# ggsave(
#   filename = "docs/img/basic/search_map_two_samples_in_one_plot.png",
#   plot = p,
#   scale = 2.5, width = 1000, height = 400, units = "px"
# )
