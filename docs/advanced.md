# Advanced features of the mobest package

`mobest::locate()` uses spatiotemporal interpolation to calculate spatial similarity probability maps between a set of search samples and an interpolated ancestry field at specific time slices. This basic and arguably most important usecase of the mobest package is documented in {doc}`A basic similarity search workflow <basic>`. `locate()` hides a lot of the complexity of mobest, though, and some of that will be documented in this section.

## Gaussian process regression on top of a linear model

An important detail of mobest is that it typically performs the spatiotemporal interpolation for an individual dependent variable in a two-step process. It starts with a simple linear model that predicts the dependent variable based on the three independent space and time variables, thus detrending the entire dataset. The more advanced, complex and computationally expensive Gaussian process regression model is then only applied to the residuals of the linear models.

The detrending step can be turned off for individual variables in `mobest::create_kernel`, by setting `on_residuals = FALSE` (see {ref}`Kernel parameter settings <types:kernel parameter settings>`).

## Spatiotemporal interpolation permutations in a model grid

Spatiotemporal interpolation is the first main operation mobest undertakes. The actual interpolation is performed in an internal function `mobest:::interpolate()`, which has a minimal interface and does not have to be called directly by the user. This is not least, because the main workflow in `mobest::locate()` (see {doc}`A basic similarity search workflow <basic>`) generally requires multiple interpolation runs. mobest therefore offers a slightly more complex, two step interpolation workflow, which consists of the creation of a list of models and then subsequently running each element in this list to construct different ancestry fields. 

`mobest::create_model_grid` creates an object of class `mobest_modelgrid` which holds all permutations of the field-defining input objects. Each row equals one complete model definition with all parameters and input data fully defined. Here is this function is called with the necessary input data constructors. Note that unlike for `locate()` we here require the `*_multi` constructors to express iterations of all settings. This mechanism is explained below in {ref}`Similarity search with permutations <advanced:similarity search with permutations>`, so we omit the details here with `...`.

```r
library(magrittr)
model_grid <- mobest::create_model_grid(
  independent     = mobest::create_spatpos_multi(...),
  dependent       = mobest::create_obs_multi(...),
  kernel          = mobest::create_kernset_multi(...),
  prediction_grid = mobest::create_spatpos_multi(...)
)
```

This model grid is only a specification of the models we want to run, so the interpolated fields we want to create. It features the following columns, some of which are list columns including the concrete data to run inform the desired interpolation.

|Column               |Description |
|:--------------------|:-----------|
|independent_table_id |Identifier of the spatiotemporal position permutation|
|dependent_setting_id |Identifier of the dependent variable space position permutation|
|dependent_var_id     |Identifier of the dependent variable|
|kernel_setting_id    |Identifier of the kernel setting permutation|
|pred_grid_id         |Identifier of the spatiotemporal prediction grid|
|independent_table    |List column: tibble with spatiotemporal positions<br>(`mobest_spatiotemporalpositions`)|
|dependent_var        |List column: numerical vector with dependent variable values|
|kernel_setting       |List column: list with kernel parameters (`mobest_kernel`)|
|pred_grid            |List column: tibble with spatiotemporal positions<br>(`mobest_spatiotemporalpositions`)|


Another function, `mobest::run_model_grid`, then takes this `model_grid` object and runs each interpolation.

```r
interpol_grid <- mobest::run_model_grid(model_grid, quiet = T)
```

It returns an unnested table of type `mobest_interpolgrid`, where each row documents the result of the interpolation for one prediction grid point and each parameter setting. It includes the following columns/variables: `independent_table_id`, `dependent_setting_id`, `dependent_var_id`, `kernel_setting_id`, `pred_grid_id`, `dsx`, `dsy`, `dt`, `g`, `id`, `x`, `y`, `z`, `mean`, `sd`. These are identical to what we already know from the output of `mobest::locate()` (see {ref}`The mobest_locateoverview table <basic:the mobest_locateoverview table>`) or `mobest::crossvalidate()` (see {ref}`A basic crossvalidation setup <estimation:a basic crossvalidation setup>`). This is simply, because these functions call `create_model_grid()` and `run_model_grid()` and build their output on top of the `mobest_interpolgrid` structure.

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
