# Advanced features of the mobest package

`mobest::locate()` uses spatiotemporal interpolation to calculate spatial similarity probability maps between a set of search samples and an interpolated ancestry field at specific time slices. This basic and arguably most important usecase of the mobest package is documented in {doc}`A basic similarity search workflow <basic>`. `locate()` hides a lot of the complexity of mobest, though, and some of that will be documented in this section.

## Gaussian process regression on top of a linear model

...

## Spatiotemporal interpolation permutations in a model grid

The spatiotemporal interpolation workflow consists of the creation of a list of models and then subsequently running each element in this list to construct different ancestry fields. The actual interpolation is done in a function `mobest:::interpolate`, which has a minimal interface and is therefore kept internal.

`mobest::create_model_grid` creates an object of class `mobest_modelgrid` which holds all permutations of the field-defining input objects. Each row equals one complete model definition with all parameters and input data fully defined.

```r
library(magrittr)
model_grid <- mobest::create_model_grid(
  independent = mobest::create_spatpos_multi(
    dating_1 = positions %>% dplyr::mutate(z = z + sample(-100:100, 100)),
    dating_2 = positions %>% dplyr::mutate(z = z + sample(-100:100, 100))
  ),
  dependent = mobest::create_obs_multi(
    obs1 = observations,
    obs2 = observations
  ),
  kernel = mobest::create_kernset_multi(
    kernel_1 = mobest::create_kernset(
      ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
      ac2 = mobest::create_kernel(1000000, 1000000, 200, 0.1)
    ),
    kernel_2 = mobest::create_kernset(
      ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
      ac2 = mobest::create_kernel(1000000, 1000000, 250, 0.1)
    )
  ),
  prediction_grid = mobest::create_spatpos_multi(
    pred_grid_1 = expand.grid(
      x = seq(100000, 1000000, 100000), 
      y = seq(100000, 1000000, 100000),
      z = seq(-5500, -3000, 500)
    ) %>% { mobest::create_spatpos(id = 1:nrow(.), x = .$x, y = .$y, z = .$z) },
    pred_grid_2 = expand.grid(
      x = seq(100000, 1000000, 100000), 
      y = seq(100000, 1000000, 100000),
      z = seq(-5500, -3000, 500)
    ) %>% { mobest::create_spatpos(id = 1:nrow(.), x = .$x, y = .$y, z = .$z) }
  )
)
```

`mobest::run_model_grid` runs each model and returns an unnested table of interpolation results for each prediction grid point and each model parameter setting.

```r
interpol_grid <- mobest::run_model_grid(model_grid, quiet = T)
```

## Similarity search with permutations

It requires the necessary reference sample input to perform the interpolation, which internally employs `mobest::create_model_grid` and `mobest::run_model_grid`. The search then yields a similarity probability value for each grid cell and for each search sample in an object of class `mobest_locateoverview`. These probabilities are normalized for each search sample and grid (with the default `normalize = TRUE`).

`mobest::locate` is actually just a special, simplified interface to `mobest::locate_multi`, which adds another level of complexity. It allows multiple input values for `independent`, `dependent`, `kernel`, `search_independent` and `search_dependent` and the result will therefore consider all permutations of these input settings (`independent` and `search_independent` as well as `dependent` and `search_dependent` have to be congruent, though).

```r
locate_overview <- mobest::locate_multi(
  independent = mobest::create_spatpos_multi(
    dating1 = positions %>% dplyr::mutate(z = z + sample(-100:100, 100)),
    dating2 = positions %>% dplyr::mutate(z = z + sample(-100:100, 100))
  ),
  dependent = mobest::create_obs_multi(
    obs1 = observations, obs2 = observations
  ),
  kernel = mobest::create_kernset_multi(
    kernel_1 = mobest::create_kernset(
      ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
      ac2 = mobest::create_kernel(1000000, 1000000, 200, 0.1)
    ),
    kernel_2 = mobest::create_kernset(
      ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
      ac2 = mobest::create_kernel(1000000, 1000000, 250, 0.1)
    )
  ),
  search_independent = mobest::create_spatpos_multi(
    dating1 = positions[1:4,], dating2 = positions[1:4,]
  ),
  search_dependent = mobest::create_obs_multi(
    obs1 = observations[1:4,], obs2 = observations[1:4,]
  ),
  search_space_grid = expand.grid(
      x = seq(100000, 1000000, 100000), 
      y = seq(100000, 1000000, 100000)
    ) %>% { mobest::create_geopos(id = 1:nrow(.), x = .$x, y = .$y) },
  search_time = 0,
  search_time_mode = "relative",
  quiet = T
)

locate_product <- mobest::multiply_dependent_probabilities(locate_overview)
```

`mobest::locate_multi` produces many probability grids for each sample. Even after `mobest::multiply_dependent_probabilities` merges the per-ancestry component iterations, that still leaves many parameter permutations. `mobest::fold_probabilities_per_group` is a convenient function to combine these to a single, merged probability grid of class `mobest_locatefold`. The folding operation can be set in the argument `folding_operation`, where the default is a simple sum. Again, in the default setting, the output probabilities are normalized per permutation.

```r
mobest::fold_probabilities_per_group(locate_product)
```

`fold_probabilities_per_group` also allows to maintain the the permutation groups, in case a full summary is not desired:

```r
mobest::fold_probabilities_per_group(locate_product, dependent_setting_id, kernel_setting_id)
```
