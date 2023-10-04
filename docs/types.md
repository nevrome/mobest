# Data types in the mobest R package

The following guide briefly lists the main `mobest` data types with their constructors loosely in the order one would usually call them.

## Basic data types

`mobest` employs a number of basic [S3 data types](http://adv-r.had.co.nz/S3.html) to formalize the input to almost all of its functions. Most of then are tabular and inherit from `tibble::tibble()`. The constructors check certain properties to insure input correctness.

### Spatial coordinates

`mobest::create_geopos` creates an object of class `mobest_spatialpositions` which is a `tibble` that represents spatial positions. Spatial positions in `mobest` are always 2-dimensional coordinates in a Cartesian space. For real world coordinates that means, that they have to be transformed to a projected coordinate system (e.g. with `sf::st_transform`): `mobest` can not be used with longitude and latitude coordinates (also see See {ref}`Projecting the spatial data <basic:projecting the spatial data>`).

```r
mobest::create_geopos(
  id = 1:100,
  x = c(sample(100000:700000, 50), sample(300000:1000000, 50)),
  y = c(sample(100000:700000, 50), sample(300000:1000000, 50))
)
```

For the interpolation fields we often want regular, spatial grids covering a specific spatial area. These can be constructed with `mobest::create_prediction_grid`, which takes an object of class `sf` with polygons in a projected coordinate system. It also yields an object of class `mobest_spatialpositions`.

Here is an example for the landmass of Europe, covered in a 250km grid:

```r
rnaturalearthdata::countries50 %>%
  sf::st_as_sf() %>%
  sf::st_make_valid() %>%
  sf::st_crop(xmin = -10.8, ymin = 33.6, xmax = 34.5, ymax = 61.3) %>%
  sf::st_transform(3857) %>%
  mobest::create_prediction_grid(250000) %>%
  ggplot() +
    geom_raster(aes(x, y)) +
    geom_text(aes(x,y,label = id), colour = "white", size = 2.5) +
    coord_fixed()
```

### Spatiotemporal coordinates

`mobest_spatialpositions` can be transformed to `mobest_spatiotemporalpositions` with `mobest::geopos_to_spatpos`. This function calculates the permutations of all spatial positions with a new, implicitly temporal dimension `z`. `mobest_spatiotemporalpositions` is also derived from `tibble`.

```r
mobest::geopos_to_spatpos(
  mobest::create_geopos(
    id = 1:100,
    x = c(sample(100000:700000, 50), sample(300000:1000000, 50)),
    y = c(sample(100000:700000, 50), sample(300000:1000000, 50))
  ),
  c(-5000, -4900, -4800)
)
```

`mobest::create_spatpos` directly creates `mobest_spatiotemporalpositions` objects to represent spatiotemporal positions.

```r
positions <- mobest::create_spatpos(
  id = 1:100,
  x = c(sample(100000:700000, 50), sample(300000:1000000, 50)),
  y = c(sample(100000:700000, 50), sample(300000:1000000, 50)),
  z = c(sample(-5000:-3500, 50), sample(-4500:-3000, 50))
)
```

### Genetic coordinates

`mobest::create_obs` creates an object `mobest_observations`, which is a `tibble` with genetic coordinates. Genetic coordinates can be any simple numeric measure of ancestry, for example the position of the samples in PCA space.

```r
observations <- mobest::create_obs(
  ac1 = c(runif(50, 0, 0.6), runif(50, 0.4, 1)), # "ac" here for "ancestry component", e.g. PCA coordinate 1
  ac2 = c(runif(50, 0, 0.3), runif(50, 0.5, 1)) # e.g. PCA coordinate 2
)
```

Names and number of the components are freely selectable, so instead of `ac1` + `ac2` as in the example here, one could, for example, also have `PC1` + `PC2` + `PC3`, or `MDS1` + `MDS2`.

### Kernel parameter settings

Gaussian process regression requires a parametrized covariance function: a "kernel". One `mobest_kernel` can be constructed with `mobest::create_kernel`. `mobest_kernel` only represents one specific kernel, though, for one specific ancestry component (e.g. `ac1`). Given that an analysis typically involves multiple genetic coordinates `mobest::create_kernset` provides a wrapper to bundle multiple kernels directly in an object of class `mobest_kernelsetting`.

```r
kernset <- mobest::create_kernset(
  ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
  ac2 = mobest::create_kernel(1000000, 1000000, 200, 0.1)
)
```

If a function requires both input of type `mobest_observations` and `mobest_kernelsetting`, then the names of the individual ancestry components must be identical, i.e. fit to each other.

## Variability and permutations

When working with real data we often need to explore permutations of data or account for uncertainty by sampling from distributions (e.g. uncertain dating). To represent that, `mobest` provides wrapper classes and constructors with a `*_multi` suffix, to bundle multiple individual elements in a list class. Some of the core functions provide interfaces that automatically consider all permutations of these input lists.

Available are:

- `mobest_spatialpositions_multi` (`mobest::create_geopos_multi`)
- `mobest_spatiotemporalpositions_multi` (`mobest::create_spatpos_multi`)
- `mobest_observations_multi` (`mobest::create_obs_multi`)
- `mobest_kernelsetting_multi` (`mobest::create_kernset_multi`)

And here is an example how they can be filled with named arguments:

```r
multiple_kernel_settings <- mobest::create_kernset_multi(
  kernel_1 = mobest::create_kernset(
    ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
    ac2 = mobest::create_kernel(1000000, 1000000, 200, 0.1)
  ),
  kernel_2 = mobest::create_kernset(
    ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
    ac2 = mobest::create_kernel(1000000, 1000000, 250, 0.1)
  )
)
```
