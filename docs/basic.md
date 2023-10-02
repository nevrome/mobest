# A basic similarity search workflow

This section explains the setup for a basic ancestry similarity search with mobest in R. 

We use a simplified version of the data and code used for the publication that introduced mobest, {cite:p}`Schmid2023`. The script explained in the following sections as well as the data required for it can be downloaded in its entirety here:

- [simple_similarity_search.R](data/simple_similarity_search.R)
- [epsg3035.?]()
- [extended_area.?]()
- [input.csv]()

## Preparing the computational environment

For this script we use various packages beyond base R, among which the following ones are required:

- `readr` for loading .csv input data
- `magrittr` for the pipe operator `%>%`
- `sf` for loading and manipulating spatial data
- `rnaturalearth` for downloading geodata
- `ggplot2` to visualize intermediate and final results
- `dplyr` for data manipulation of `data.frame`s
- `mobest` (obviously)

`readr`, `magrittr`, `ggplot2` and `dplyr` are all in the [tidyverse](https://www.tidyverse.org) and can be installed in one go with `install.packages("tidyverse")` on the R console.

For the installation of `sf` and `mobest` please see the instructions in the [Install section](install.md).

We will generally call functions explicitly with their namespace using `::` (so e.g. `readr::read_csv()`). The only exceptions are `magrittr` and `ggplot2`, because we will use their functions so often that it becomes tedious to type them out. Instead we load them at the beginning.

```r
library(magrittr)
library(ggplot2)
```

## Preparing the input data

### Generating the the spatial prediction grid

mobest's similarity search is typically run for a regular grid of spatial positions in the area of interest. It provides a function (`mobest::create_prediction_grid()`) to create such a grid, given a specification of the desired area. This area is typically (based on how we imagine mobest to be used) the land area in a certain part of planet Earth.

#### Defining the research area

In a first step we therefore have to define the research area for our analysis as a polygon in space. One way of doing this is to provide a list of latitude and longitude coordinates (extracted e.g. from Google Maps). The following code defines a simple research area covering large parts of Western Eurasia.

```r
research_area_4326 <- sf::st_polygon(
  list(
    cbind(
      c(35.91,11.73,-11.74,-15.47,37.06,49.26,49.56,35.91), # longitudes
      c(25.61,28.94, 31.77, 62.73,65.67,44.56,28.55,25.61)  # latitudes
    )
  )
) %>% sf::st_sfc(crs = 4326)
```

Spatial coordinates require a coordinate references system (CRS). For lat-lon coordinates we typically use [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System#WGS84) with the [EPSG code](https://en.wikipedia.org/wiki/EPSG_Geodetic_Parameter_Dataset) 4326. `st_polygon()` creates a simple polygon as a clockwise arrangement of individual coordinates and `st_sfc()` properly defines this polygon as a geographic area on Earth. A simple way to interactively inspect this polygon on a world map in R is provided by the mapview package: `mapview::mapview(research_area_4326)`.

```{figure} img/basic/mapview_research_area.png
The defined research area plotted on top of a map.
```

With the research area properly defined we can move to the next challenge and extract the land area in the research area. For that we have to obtain a dataset with polygons that trace the coastlines for all around the world. The [naturalearthdata](https://www.naturalearthdata.com) project provides open worldwide geodata in different resolutions and in easy to use data formats. The `rnaturalearth` package makes it easy to download this data right into `sf` objects in R.

```r
worldwide_land_outline_4326 <- rnaturalearth::ne_download(
  scale = 50, type = 'land', category = 'physical',
  returnclass = "sf"
)
```

We can then crop the land outline to the research area to obtain the land area we are interested in.

```r
research_land_outline_4326 <- sf::st_intersection(
  worldwide_land_outline_4326,
  research_area_4326
)
```

We can plot the resulting spatial multi-polygon with ggplot2.

```r
ggplot() + geom_sf(data = research_land_outline_4326)
```

```{figure} img/basic/research_area_land_outline_4326.png
The research area land polygon.
```

#### Projecting the spatial data

At this point we run into a specific issue of mobest: It requires its "independent" spatial and temporal coordinates to be coordinates in a [cartesian system](https://en.wikipedia.org/wiki/Cartesian_coordinate_system) describing [Euclidean space](https://en.wikipedia.org/wiki/Euclidean_space). For the spatial coordinates that means we can not work with latitude and longitude coordinates on a sphere, but have to apply [map projection](https://en.wikipedia.org/wiki/Map_projection) to represent the curved, two dimensional surface of our planet on a simple plane.

The question how exactly this should be done and which CRS to choose depends on the position, size and shape of your research area. Each map projection algorithm has different properties regarding whether they manage to preserve or distort size, shape, distances and directions of areas and lines compared to the actual properties on Earth. Generally the larger the research area, the bigger the distortion of these properties becomes and for mobest we ideally want to represent all them accurately. mobest is therefore unfit for origin search on a global scale, but can usually be well applied for individual countries with the projections recommended by their cartographic agencies. For an intermedite, continental scale, as in this example, we have to choose our CRS wisely. 

We decided to follow the recommendation of {cite:p}`Annoni2003`.

> The Workshop recommends that the European Commission:    
> Uses for statistical analysis and display a ETRS89 Lambert Azimuthal Equal Area coordinate reference system of 2001 [ETRS -LAEA11 ], that is specified by ETRS89 as datum and the Lambert Azimuthal Equal Area map projection.    
> ...

This setting is documented in the EPSG code [3035](https://epsg.io/3035). Our decision comes at the price of increased inaccuracy especially in the North- and South-East of the research area where we get very far away from the center at 52° latitude and 10° longitude (see {cite:p}`Tsoulos2003` p.53 for a visualization of the deformation).

To transform the the land outline in the research area from EPSG:4326 to EPSG:3035 we can apply `sf::st_transform()`.

```r
research_land_outline_3035 <- research_land_outline_4326 %>% sf::st_transform(3035)
```

Note how the change in the coordinate system affects the map plot.

```r
ggplot() + geom_sf(data = research_land_outline_3035)
```

```{figure} img/basic/research_area_land_outline_3035.png
The research area land polygon now transformed to EPSG:3035.
```

#### Creating the prediction grid

```r
spatial_pred_grid <- mobest::create_prediction_grid(
  research_land_outline_3035,
  spatial_cell_size = 50000
)
```

```r
ggplot() +
  geom_point(data = spatial_pred_grid, mapping = aes(x, y), color = "red")
```

### Reading the input samples

## Running mobest's interpolation and search function

### Building the input data structures

### Calling `mobest::locate`

## Inspecting the computed results

<!--

## Artifical example

Here is a simple, artificial example how 2. can be used:

```r


# a function to calculate the similarity probability for one particular sample
locate_simple <- mobest::locate(
  # spatiotemporal coordinates of the reference samples informing the ancestry field
  independent = mobest::create_spatpos(
    id = 1:100,
    x = c(sample(100000:700000, 50), sample(300000:1000000, 50)), # space x
    y = c(sample(100000:700000, 50), sample(300000:1000000, 50)), # space y
    z = c(sample(-5000:-3500, 50), sample(-4500:-3000, 50))       # time
  ),
  # genetic coordinates of the reference samples
  dependent = mobest::create_obs(
    ac1 = c(runif(50, 0, 0.6), runif(50, 0.4, 1)), # PCA coordinate 1
    ac2 = c(runif(50, 0, 0.3), runif(50, 0.5, 1))  # PCA coordinate 2
  ),
  # field properties for each ancestry component
  kernel = mobest::create_kernset(
    ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
    ac2 = mobest::create_kernel(1000000, 1000000, 200, 0.1)
  ),
  # spatiotemporal coordinates of the sample of interest
  search_independent = mobest::create_spatpos(
    id = 1,
    x = sample(100000:1000000, 1), # space x
    y = sample(100000:1000000, 1), # space y
    z = sample(-5000:-3000, 1)     # time
  ),
  # genetic coordinates of the sample of interest
  search_dependent = mobest::create_obs(
    ac1 = runif(1, 0, 0.6), # PCA coordinate 1
    ac2 = runif(1, 0, 0.5)  # PCA coordinate 2
  ),
  # spatial search grid: Where to search
  search_space_grid = expand.grid(
      x = seq(100000, 1000000, 100000), 
      y = seq(100000, 1000000, 100000)
    ) %>% { mobest::create_geopos(id = 1:nrow(.), x = .$x, y = .$y) },
  # search time: When to search
  search_time = -500,
  quiet = T
)

# multiply probabilities for PCA coordinate 1 and PCA coordinate 2
locate_product <- mobest::multiply_dependent_probabilities(locate_simple)

# plot the resulting probability surface
library(ggplot2)
locate_product %>% ggplot() +
  geom_raster(mapping = aes(x = field_x, y = field_y, fill = probability)) +
  geom_point(mapping = aes(x = search_x, y = search_y), colour = "red") +
  coord_fixed() +
  ggtitle(paste0(
    "t for sample of interest = ", unique(locate_product$search_z), "\n",
    "t field time slice = ", unique(locate_product$field_z)
  ))
```

## Origin search

`mobest::locate` uses the spatiotemporal interpolation to calculate a similarity probability between a set of "search" samples and an interpolation field. It requires the necessary reference sample input to perform the interpolation, which internally employs `mobest::create_model_grid` and `mobest::run_model_grid`. The search then yields a similarity probability value for each grid cell and for each search sample in an object of class `mobest_locateoverview`. These probabilities are normalized for each search sample and grid (with the default `normalize = TRUE`).

```r
locate_simple <- mobest::locate(
  independent = positions,
  dependent = observations,
  kernel = kernset,
  search_independent = positions[1:4,],
  search_dependent = observations[1:4,],
  search_space_grid = expand.grid(
      x = seq(100000, 1000000, 100000), 
      y = seq(100000, 1000000, 100000)
    ) %>% { mobest::create_geopos(id = 1:nrow(.), x = .$x, y = .$y) },
  search_time = c(0,-100),
  quiet = T
)
```

The spatiotemporal probability grids `locate` returns are calculated are per ancestry component (as put in via `dependent`/`search_dependent`). To multiply the ancestry component grids, `mobest` provides `mobest::multiply_dependent_probabilities`, which yields an object of class `mobest_locateproduct`. The output probabilities are normalized per permutation.

```r
mobest::multiply_dependent_probabilities(locate_simple)
```

-->
