library(magrittr)
library(ggplot2)

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

# ggplot() + geom_sf(data = research_land_outline_4326)

research_land_outline_3035 <- research_land_outline_4326 %>% sf::st_transform(crs = 3035)

# ggplot() + geom_sf(data = research_land_outline_3035)

spatial_pred_grid <- mobest::create_prediction_grid(
  research_land_outline_3035,
  spatial_cell_size = 50000
)

# ggplot() +
#   geom_sf(data = research_land_outline_3035) +
#   geom_point(
#     data = spatial_pred_grid,
#     mapping = aes(x, y),
#     color = "red",
#     size = 0.25
#   )

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

# ggplot() +
#   geom_sf(data = research_land_outline_3035) +
#   geom_point(
#     data = samples_projected,
#     mapping = aes(x, y),
#     color = "darkgreen",
#     size = 0.25
#   )

search_samples <- samples_projected %>%
  dplyr::filter(
    Sample_ID == "Stuttgart_published.DG"
  )

search_result <- mobest::locate(
  independent = mobest::create_spatpos(
    id = samples_projected$Sample_ID,
    x  = samples_projected$x,
    y  = samples_projected$y,
    z  = samples_projected$Date_BC_AD_Median
  ),
  dependent = mobest::create_obs(
    C1 = samples_projected$MDS_C1,
    C2 = samples_projected$MDS_C2
  ),
  kernel = mobest::create_kernset(
    C1 = mobest::create_kernel(
      dsx = 800 * 1000, dsy = 800 * 1000, dt = 800,
      g = 0.1
    ),
    C2 = mobest::create_kernel(
      dsx = 800 * 1000, dsy = 800 * 1000, dt = 800,
      g = 0.1
    )
  ),
  search_independent = mobest::create_spatpos(
    id = search_samples$Sample_ID,
    x  = search_samples$x,
    y  = search_samples$y,
    z  = search_samples$Date_BC_AD_Median
  ),
  search_dependent = mobest::create_obs(
    C1 = search_samples$MDS_C1,
    C2 = search_samples$MDS_C2
  ),
  search_space_grid = spatial_pred_grid,
  search_time = -6500,
  search_time_mode = "absolute"
)
