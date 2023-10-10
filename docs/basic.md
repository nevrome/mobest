# A basic similarity search workflow

This section explains the setup for a basic ancestry similarity search with mobest in R. 

For didactic purposes we use a simplified version of the data and code used for the publication that introduced mobest: {cite:p}`Schmid2023`. This is a fairly generic setup you can adjust to the needs of other and more specific projects.

The script explained in the following sections as well as the data required for it can be downloaded in its entirety here:

- [simple_similarity_search.R](data/simple_similarity_search.R)
- [samples_basic.csv](data/samples_basic.csv)

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

At this point we run into a specific issue of mobest: It requires its "independent" spatial and temporal coordinates to be coordinates in a [Cartesian system](https://en.wikipedia.org/wiki/Cartesian_coordinate_system) describing [Euclidean space](https://en.wikipedia.org/wiki/Euclidean_space). For the spatial coordinates that means we can not work with latitude and longitude coordinates on a sphere, but have to apply [map projection](https://en.wikipedia.org/wiki/Map_projection) to represent the curved, two dimensional surface of our planet on a simple plane.

The question how exactly this should be done and which CRS to choose depends on the position, size and shape of your research area. Each map projection algorithm has different properties regarding whether they manage to preserve or distort size, shape, distances and directions of areas and lines compared to the actual properties on Earth. Generally the larger the research area, the bigger the distortion of these properties becomes and for mobest we ideally want to represent all them accurately. mobest is therefore unfit for origin search on a global scale, but can usually be well applied for individual countries with the projections recommended by their cartographic agencies. For an intermediate, continental scale, as in this example, we have to choose our CRS wisely. 

We decided to follow the recommendation of {cite:p}`Annoni2003`.

> The Workshop recommends that the European Commission:    
> Uses for statistical analysis and display a ETRS89 Lambert Azimuthal Equal Area coordinate reference system of 2001 [ETRS -LAEA11 ], that is specified by ETRS89 as datum and the Lambert Azimuthal Equal Area map projection.    
> ...

This setting is documented in the EPSG code [3035](https://epsg.io/3035). Our decision comes at the price of increased inaccuracy especially in the North- and South-East of the research area where we get very far away from the center at 52° latitude and 10° longitude (see {cite:p}`Tsoulos2003` p.53 for a visualization of the deformation).

To transform the the land outline in the research area from EPSG:4326 to EPSG:3035 we can apply `sf::st_transform()`.

```r
research_land_outline_3035 <- research_land_outline_4326 %>% sf::st_transform(crs = 3035)
```

Note how the change in the coordinate system affects the map plot.

```r
ggplot() + geom_sf(data = research_land_outline_3035)
```

```{figure} img/basic/research_area_land_outline_3035.png
The research area land polygon now transformed to EPSG:3035.
```

#### Creating the prediction grid

To finally create the prediction grid we can use `mobest::create_prediction_grid()`. It takes the land outline polygon and overlays its bounding box with a regular grid (using `sf::st_make_grid()`), where each cell has the size corresponding to the `spatial_cell_size` parameter. It then determines the centers of each grid cell and crops the resulting, regular point cloud with the land area. Note that `spatial_cell_size` uses the unit of the CRS, so in our case for EPSG:3035 meters. That means a value of 50000 translates to one point every 50km. The total number of resulting spatial prediction positions is 4738 in this example.

```r
spatial_pred_grid <- mobest::create_prediction_grid(
  research_land_outline_3035,
  spatial_cell_size = 50000
)
```

`create_prediction_grid` returns an object of class `mobest_spatialpositions`, which is derived from `tibble::tibble`. That means we can print it on the R console and it will behave as a tibble. It will also work seamlessly as an input for ggplot2, which we can now use to visualize the point cloud.

```r
ggplot() +
  geom_sf(data = research_land_outline_3035) +
  geom_point(
    data = spatial_pred_grid,
    mapping = aes(x, y),
    color = "red",
    size = 0.25
  )
```

```{figure} img/basic/spatial_prediction_grid.png
The spatial prediction grid points plotted on top of the land area.
```

### Reading the input samples

mobest requires a set of data points, observations, to inform the ancestry field interpolation. For each observation the position in space, time and a dependent variable space (e.g. the coordinates in a PCA analysis) must be known. This information must be provided in a specific format. A typical workflow would involve preparing this information in a .xlsx or (better) .csv table, which could then be read into R.

For this tutorial we rely on the data used and published in {cite:p}`Schmid2023`. The following, hidden section includes the code to prepare the sample table we need here directly from the supplementarty tables published with the paper.

<details>
<summary>Code to prepare the input data table.</summary>

```r
# download .zip archives with tables from https://doi.org/10.17605/OSF.IO/6UWM5
utils::download.file(
  url = "https://osf.io/download/kej4s/",
  destfile = "docs/data/pnas_tables.zip"
)
# extract the relevant tables
utils::unzip(
  "docs/data/pnas_tables.zip",
  files = c("Dataset_S1.csv", "Dataset_S2.csv"),
  exdir = "docs/data/"
)
# read data files
samples_context_raw <- readr::read_csv("docs/data/Dataset_S1.csv")
samples_genetic_space_raw <- readr::read_csv("docs/data/Dataset_S2.csv")
# join them by sample name
samples_raw <- dplyr::left_join(
  samples_context_raw,
  samples_genetic_space_raw,
  by = "Sample_ID"
)
# create different useful subsets of this table
## most basic selection of variables
samples_basic <- samples_raw %>%
  dplyr::select(
    Sample_ID,
    Latitude, Longitude,
    Date_BC_AD_Median,
    MDS_C1 = C1_mds_u, MDS_C2 = C2_mds_u
  )
readr::write_csv(samples_basic, file = "docs/data/samples_basic.csv")
```

</details>

You do not have to run this and can instead download the example table `samples_basic.csv` from [here](data/samples_basic.csv).

When you have download the input data file you can load it into a [`tibble`](https://tibble.tidyverse.org) in R.

```r
samples_basic <- readr::read_csv("docs/data/samples_basic.csv")
# you have to replace "data/docs/" with the path to your file
```

`samples_basic` is a tibble that includes the following columns/variables:

| Variable          | Type | Description                                                                                                              |
|-------------------|------|--------------------------------------------------------------------------------------------------------------------------|
| Sample_ID         | chr  | A sample identifier                                                                                                      |
| Latitude          | dbl  | The latitude coordinate where this sample was recovered                                                                  |
| Longitude         | dbl  | The longitude coordinate                                                                                                 |
| Date_BC_AD_Median | int  | The median age of this sample in years                                                                                   |
| MDS_C1            | dbl  | The coordinate of this sample on dimension 1 of an MDS analysis.<br>See the paper for more details on how this was obtained |
| MDS_C2            | dbl  | The coordinate of this sample on MDS dimension 2                                                                         |

These variables are a minimum for a meaningful mobest run and must be known for all samples. Samples that are missing information in any of these columns have to excluded from the input data table.

Before we move on, we have to apply one more change to the sample table: Just as for the research area (see {ref}`Projecting the spatial data <basic:projecting the spatial data>`) above) we have to transform the coordinates from longitude and latitude coordinates to a projected system, specifically the same we selected above. To do this we can construct an `sf` object from the sample `tibble`, apply `sf::st_transform` and then transform this result back to a `tibble` with the x and y coordinates of EPSG:3035 in extra columns. This last step makes the code a bit awkward.


```r
samples_projected <- samples_basic %>%
  # make the tibble an sf object
  sf::st_as_sf(
    coords = c("Longitude", "Latitude"),
    crs = 4326
  ) %>%
  # transform the coordinates
  sf::st_transform(crs = 3035) %>% 
  # reshape the sf object back into a simple tibble
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],
    y = sf::st_coordinates(.)[,2]
  ) %>%
  sf::st_drop_geometry()
```

With the coordinates in the same reference system as the land outline we prepared above we can combine both in a figure:

```
ggplot() +
  geom_sf(data = research_land_outline_3035) +
  geom_point(
    data = samples_projected,
    mapping = aes(x, y),
    color = "darkgreen",
    size = 0.25
  )
```

```{figure} img/basic/samples_map.png
The spatial distribution of the informative samples.
```

A number of samples are outside of the area we actually want to predict here. That is no problem. They will inform the field in the north-eastern fringes of the area of interest and do no harm. It is much more problematic that some areas of our prediction grid are severely under-sampled. That is something we have to keep in mind for later when we interpret the results of the similarity search.

## Specifying the search sample

mobest's similarity search always takes the perspective of an individual sample for which we want to determine similarity probabilities for a spatial prediction grid at a specific time. For this sample, the "search sample" we require the same information as for the input samples: The position in space, time and the dependent variable space (e.g. PCA or MDS space).

Technically this is only a requirement of the mobest interface. Conceptually such a similarity search only really requires the dependent variable space position of interest. The added benefit of having all information there is the relative time setting (see below) and a very comprehensive output table for the typical use-case.

In this example we will use one specific sample with a pretty well studied ancestry history: The sample named `Stuttgart` published in {cite}`Lazaridis2014`. We can select it as a subset of our sample table:

```r
search_samples <- samples_projected %>%
  dplyr::filter(
    Sample_ID == "Stuttgart_published.DG"
  )
```

With this setup the search sample itself will be part of the samples used to inform the ancestry field interpolation. This is no problem - the search sample is a known data point in space and time that can very well be employed to build a better model of the past ancestry distribution. There may be research questions for which this might not be desired, though. Then it can just as well be excluded from the `samples_projected` table.

## Running mobest's interpolation and search function

With the input data, both the spatial prediction grid and the samples to inform the ancestry field interpolation, prepared and ready, we can now run `mobest::locate`. For that we first have to split and transform the input data into the required data structures.

### Building the input data for the interpolation

Here is how the interface of `mobest::locate` looks:

```r
mobest::locate(
  # spatiotemporal coordinates of the reference samples informing the ancestry field
  independent = ...,
  # genetic coordinates of the reference samples
  dependent   = ...,
  # ---
  # interpolation settings for each ancestry component
  kernel = ...,
  # ---
  # spatiotemporal coordinates of the sample of interest
  search_independent = ...,
  # genetic coordinates of the sample of interest
  search_dependent   = ...,
  # ---
  # spatial search grid: Where to search
  search_space_grid  = ...,
  # search time: When to search
  search_time      = ...,
  search_time_mode = ...,
  # ---
  # should the result be normalized
  normalize = ...
)
```

Each of these arguments requires specific input.

#### Independent and dependent positions

The `locest()` arguments `independent` and `dependent` take the spatiotemporal and "genetic" (e.g. MDS/PCA) positions of the interpolation-informing samples. The terms *independent* and *dependent* allude to the notion and terminology of a statistical model, where positions in dependent, genetic space are predicted based on positions in independent, spatiotemporal space.

Spatiotemporal positions are encoded in `mobest` with a custom data type: {ref}`mobest_spatiotemporalpositions <types:spatiotemporal coordinates>`. For the `independent` argument of `locest()` we have to construct an object of this type with `mobest::create_spatpos()` to represent the positions of the input samples in `samples_projected`.

```r
ind <- mobest::create_spatpos(
  id = samples_projected$Sample_ID,
  x  = samples_projected$x,
  y  = samples_projected$y,
  z  = samples_projected$Date_BC_AD_Median
)
```

The dependent, genetic variables are also encoded in a custom, tabular type: {ref}`mobest_observations <types:genetic coordinates>` with the constructor function `mobest::create_observations`.

```r
dep <- mobest::create_obs(
  C1 = samples_projected$MDS_C1,
  C2 = samples_projected$MDS_C2
)
```

Note that you can have an arbitrary amount of these components with arbitrary names. The only condition is, that the very same names are used below for the search samples and for the kernel parameter settings of each dependent variable.

The lengths of the vectors (`samples_projected$...`) used for `create_spatpos` and `create_obs` all have to be identical. And their order has to be the same as well: although the input is distributed over two constructors they describe the same samples.

For the search sample we have to construct objects of the same type and structure.

```r
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
```

#### Kernel parameter settings

The `locest()` argument `kernel` takes an object of the class {ref}`mobest_kernelsetting <types:kernel parameter settings>`. This type encodes kernel configurations for each dependent variable, so the parameters for the Gaussian process regression interpolation that should be used for this variable. These include mostly the lengthscale parameters in space (x and y) and time, as well as the nugget parameter. In very simple terms the former specify how far in space and time an individual sample's genetic position should inform the interpolated field. The nugget, on the other hand, is an error term to model local (for observations from the same position in space and time) variability. See {cite:p}`Gramacy2020`, specifically [here](https://bookdown.org/rbg/surrogates/chap5.html#chap5hyper), for more details.

Here is a possible configuration for our example. We construct two kernel settings, one for each ancestry component, with `mobest::create_kernel()` in `mobest::create_kernset()`. 

```
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
```

Note how we scale the lengthscale parameters: `dsx` and `dsy` are set in meters (800 * 1000m = 800km) and `dt` in years (800y). `g` is dimensionless. With the setting specified here both dependent variables will be interpolated with the same, very smooth (several hundred kilometers and years in diameter) kernel.

The main question naturally arising from this, is how to set these parameters for a given dataset and research question. There are various empirical ways to find optimal values through numerical optimization. See Supplementary Text 2 of {cite:p}`Schmid2023` for the approaches we applied. One concrete workflow to estimate the nugget from the variogram and the lengthscale parameters through crossvalidation is explained in {doc}`Interpolation parameter estimation <estimation>`.

We would argue, though, that this computationally expensive workflow is not necessary for basic applications of mobest. The analysis in {cite:p}`Schmid2023` showed that Gaussian process regression returns reasonably accurate interpolation results for a large range of kernel parameter settings, as long as they reflect a plausible intuition about the mobility behaviour of human ancestry, which generally operates on a scale of hundreds of kilometers and years. mobest is primarily a visualization method and adjusting its parameters to ones liking is legitimate if the choices are communicated transparently.

#### Search positions

With input data and settings out of the way we can now specify the points in space and time where we actually want to perform the search. For these positions the GPR model is queried to return a mean and error, which are in turn used to calculate the probability density of a specific dependent variable space position, e.g. a specific coordinate on the first coordinate of an MDS analysis.

We already performed all necessary work for the `search_space_grid` argument, so the spatial positions of the prediction grid, in {ref}`Generating the the spatial prediction grid <basic:generating the the spatial prediction grid>`. We can just enter `spatial_pred_grid` here.

The search time can be specified as an integer vector of years: e.g. `search_time = c(-500, -200, 100)`. This vector gets interpreted by `mobest::locate()` in two different ways, which can be selected with the switch argument `search_time_mode`. `search_time_mode` can either be `"relative"` (which is the default!) or `absolute`.

- With the `"relative"` mode the `search_time` is interpreted as a $\Delta t$ relative to the age of the search sample(s). Negative values point to ages that are older then the sample age, so in their relative past, and positive ones to younger ages in their relative future. In this example `-500` would be interpreted as 500 years prior to the year the Stuttgart sample presumably died (so -5242-500 = -5742 BC/AD), and 100 as an age 100 years after their death (so -5242+100 = -5142 BC/AD).
- The `"absolute"` mode is more straight forward. In this case the values in `search_time` are just interpreted as absolute ages in years BC/AD.

For this example we will set the search time to an `"absolute"` value.

```r
search_time = -6500
search_time_mode = "absolute"
```

This will search at exactly one point in time; a single timeslice 6500 BC.

#### Normalization

The last relevant option of `locate()`, `normalize`, concerns the normalization of the output. mobest's search calculates probability densities for each search point. This is a dimensionless measure that is hard to compare across multiple runs with different parameter settings. If `normalize` is set to `TRUE`, then the densities for sets of spatial points that share all other parameters (including the timeslice) are rescaled to sum to one.

We assume users generally want to use mobest, specifically `locate()`, to calculate similarity probability density maps for individual samples, time slices and parameter settings. The most natural normalization for this case is to unify the scaling of these maps. This renders them comparable.

`normalize` should therefore be set to `TRUE` for basic applications. This is also encoded as the the default setting.

### Calling `mobest::locate`

In the previous sections we have thoroughly prepared the input for a first, simple run of `mobest::locate()`. We can now call the function.

```r
search_result <- mobest::locate(
  independent        = ind,
  dependent          = dep,
  kernel             = kernset,
  search_independent = search_ind,
  search_dependent   = search_dep,
  search_space_grid  = spatial_pred_grid,
  search_time        = -6500,
  search_time_mode   = "absolute"
)
```

This typically runs for a couple of seconds, uses every available processor core and returns an object `search_result`, which we will inspect below.

## Inspecting the computed results

`mobest::locate()` returns an object of class `mobest_locateoverview`. It includes the relevant information for visualization and further processing of the analysis results.

### The mobest_locateoverview table

The output data type `mobest_locateoverview` is derived from `tibble` and has a large set of columns, many not immediatelly relevant to the basic example here. This applies especially for the variables documenting the excessive permutation mechanics hidden behind the relatively simple interface of `mobest::locate()`. `locate()` is, in fact, a wrapper function for the more flexible function `mobest::locate_multi()`, which can handle permutations in various additional input parameters (see {doc}`Advanced mobest features <advanced>`).

Spelled out this means, each row of the `mobest_locateoverview` table stores the calculated interpolated mean, error and similarity probability (`field_mean`, `field_sd`, `probability`) for one permutation of the input point positions in independent and dependent variable space (`independent_table_id` and `dependent_setting_id`), one dependent variable `dependent_var_id`, one iteration of the kernel settings (`kernel_setting_id`: `dsx`, `dsy`, `dt`, `g`), one prediction grid point emerging as a combination of spatial grid and search timeslice (`pred_grid_id`: `field_id`, `field_geo_id`, `field_x`, `field_y`, `field_z`, `search_time`) and finally one search sample (`search_id`, `search_x`, `search_y`, `search_z`, `search_measured`).
 
Here is a list of the variables returned in `mobest_observations` for each of these result iterations.

|Column               |Description |
|:--------------------|:-----------|
|independent_table_id |Identifier of the spatiotemporal position permutation|
|dependent_setting_id |Identifier of the dependent variable space position permutation|
|dependent_var_id     |Identifier of the dependent variable|
|kernel_setting_id    |Identifier of the kernel setting permutation|
|pred_grid_id         |Identifier of the spatiotemporal prediction grid|
|dsx                  |Kernel lengthscale parameter on the spatial x axis|
|dsy                  |Kernel lengthscale parameter on the spatial y axis|
|dt                   |Kernel lengthscale parameter on the temporal axis|
|g                    |Kernel nugget parameter|
|field_id             |Identifier of the spatiotemporal prediction point|
|field_x              |Spatial x axis coordinate of the prediction point|
|field_y              |Spatial y axis coordinate of the prediction point|
|field_z              |Temporal coordinate (age) of the prediction point|
|field_geo_id         |Identifier of the spatial prediction point|
|field_mean           |Mean value predicted by the GPR model for the relevant dependent variable|
|field_sd             |Uncertainty predicted by the GPR model for the relevant dependent variable|
|search_id            |Identifier of the search sample|
|search_x             |Spatial x axis coordinate of the search sample|
|search_y             |Spatial y axis coordinate of the search sample|
|search_z             |Temporal coordinate (age) of the search sample|
|search_time          |Search time as provided by the user in `locate()`'s `search_time` argument|
|search_measured      |Genetic coordinate of the search sample in the dependent variable space|
|probability          |Probability density for `search_measured` given all other parameters|

As a result of the permutation of parameters, prediction grid and search points the number of rows of `mobest_locateoverview` table can be calculated as a product of the individual counts of all relevant entities. One way to quickly validate the output of `locate()` and `locate_multi()` is to calculate the number of expected results based on the input and compare it with the actual number of rows in the output. For our example this calculation is fairly simple:

We have:

- $1$ set of input point positions in independent variable space
- $1$ set of input point positions in dependent variable space
- $2$ dependent variables
- $1$ set of kernel parameter settings
- $4738$ spatial prediction grid positions
- $1$ time slice of interest
- $1$ search sample

This means we expect exactly $2 * 4738 = 9476$ rows in `search_result`, which we can confirm with `nrow(search_result)`. 

### Creating similarity probability maps for individual dependent variables

The most basic similarity probability map we can create with `search_result` is a map for just one parameter permutation, including only one dependent variable. In this case the relevant similarity probability observations are easy to obtain. We can just filter by `dependent_var_id` to only include either `C1` or `C2`.

```r
result_C1 <- search_result %>% dplyr::filter(dependent_var_id == "C1")
```

And this is then easy to plot with `geom_raster()`. We can also plot C1 and C2 together using `cowplot::plot_grid()`.

```r
p_C1 <- ggplot() +
  geom_raster(
    data = result_C1,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  coord_fixed()

# for C2
result_C2 <- search_result %>% dplyr::filter(dependent_var_id == "C2")
p_C2 <- ggplot() +
  geom_raster(
    data = result_C2,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  coord_fixed()

# arrange both plots together
cowplot::plot_grid(p_C1, p_C2, labels = c("C1", "C2"))
```

```{figure} img/basic/search_map_simple_C1_C2.png
The similarity probability search results for the sample Stuttgart for 6500 BC.
```

### Combining the information from multiple dependent variables

The results for individual dependent variables, so ancestry components like MDS or PCA dimensions, can be informative, but are usually under-powered to exclude highly improbable search results. Generally combining them improves the accuracy of the results for individual samples, and we think this is best done by multiplying the results for the different dependent variables. This way spatial areas with high similarity probability for all dependent variables are naturally up-weighted, whereas areas that are unlikely similar for some dependent variables are down-weighted.

To perform the multiplication (and the re-normalization afterwards), mobest includes a function `mobest::multiply_dependent_probabilities()`. It works on objects of type `mobest_locateoverview` and yields tabular objects of type `mobest_locateproduct`. For this transformation it is aware of the parameter permutations potentially encoded in the `mobest_locateoverview` overview table. It only combines the probabilities for dependent variables that share all other parameters. That means the number of rows in `mobest_locateproduct` will be $\frac{1}{\text{Number of dependent variables}}$ times the number of rows in the input `mobest_locateoverview` table.

If we call it for `search_result` the output will have again $9476/2=4738$ rows. 

```r
search_product <- mobest::multiply_dependent_probabilities(search_result)
```

`mobest_locateproduct` tables have a perfect subset of the columns of `mobest_locateoverview`. We can plot the combined similarity probability map with the code already applied for individual dependent variables.

```r
ggplot() +
  geom_raster(
    data = search_product,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  coord_fixed()
```

```{figure} img/basic/search_map_simple_combined.png
The combined ($\text{C1}*\text{C2}$) similarity probability search results for the sample Stuttgart for 6500 BC.
```

This figure is the main result of this write-up.

***

## Improving the similarity search map plot

This figure is not particularly beautiful, though. We can and should apply some changes to this plot to make it more readable and visually pleasing.

1. Increase the spatial resolution of the prediction grid. `locate()` takes more time to compute with this change, but for individual samples and time slices it is generally very much affordable to go to higher resolutions. Here we go from the 50km grid above to a much finer 20km grid. The number of spatial prediction points increases from 4738 to 29583.

```r
spatial_pred_grid <- mobest::create_prediction_grid(
  research_land_outline_3035,
  spatial_cell_size = 20000
)
```

2. Change the colour scale. We chose the highly readable and colourblind-safe viridis palette.

```r
ggplot() +
  geom_raster(
    data = search_product,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  coord_fixed() +
  scale_fill_viridis_c()
```

3. Add helpful elements to the plot. We add the border of the research area, a marker for the spatial position of the search sample, a plot title and an annotation to indicate the search timeslice.

```r
research_area_3035 <- research_area_4326 %>% sf::st_transform(3035)

ggplot() +
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
    label = "6500BC",
    x = Inf, y = Inf, hjust = 1.1, vjust = 1.5
  )
```

4. Adjust plot layout details. We switch to the plot theme `theme_bw`, turn off axis labels and adjust the legend title.

```r
... +
  theme_bw() +
  theme(
    axis.title = element_blank()
  ) +
  guides(
    fill = guide_colourbar(title = "Similarity\nsearch\nprobability")
  )
```

This leaves us with the following, more refined plot:

```{figure} img/basic/search_map_neat_combined.png
A more beautiful and informative version of the result plot.
```

## Simple permutations

In the example above we performed the similarity search with `mobest::locate()` for only one parameter permutation, keeping everything constant except the two dependent variables. But as already laid out above in {ref}`The mobest_locateoverview table <basic:the mobest_locateoverview table>`, mobest can automatically consider more parameter permutations, the most basic of which are directly available in `locate()`.

Please note that all parameter permutations will be multiplied with all other permutations, causing the number of runs to grow rapidly. If you, for example, submit five time slices and five search samples, the number of runs will be $5*5=25$ times bigger than for one time slice and sample.

### Multiple search time slices

As explained in {ref}`Search positions <basic:search positions>` the `search_time` argument can take an integer vector of relative or absolute ages. That means we can run the search not just for one, but for arbitrarily many time slices at with one call to `locate`.

Here is an example with two time slices.

```r
search_result <- mobest::locate(
  independent        = ind,
  dependent          = dep,
  kernel             = kernset,
  search_independent = search_ind,
  search_dependent   = search_dep,
  search_space_grid  = spatial_pred_grid,
  search_time        = c(-6500, -5700),
  search_time_mode   = "absolute"
)
search_product <- mobest::multiply_dependent_probabilities(search_result)

ggplot() +
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
```

```{figure} img/basic/search_map_two_timeslices.png
Similarity search map plot for two different time slices.
```

### Multiple search samples

We can also select multiple search samples and prepare the input data for `mobest::locate()` 

```r
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
```

Here we set the `search_time_mode` to `"relative"` to get a different, meaningful search time for both samples.

```r
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

ggplot() +
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
```

```{figure} img/basic/search_map_two_samples.png
Similarity search map plot for two different search samples.
```