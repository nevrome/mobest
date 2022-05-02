
<!-- Rmd -> md -->

# mobest

This R package provides types and functions for spatiotemporal
interpolation of human genetic ancestry components and a derived measure
for **mob**ility **est**imation. The workflow in version 1.0 was
specifically developed to support this research compendium:
<https://github.com/nevrome/mobest.analysis.2022>.

0.  `mobest` assumes you have a set of genetic samples with spatial (two
    coordinates in a projected reference system) and temporal positions
    (years BC/AD) for which you calculated a derived, numeric measure of
    genetic ancestry (e.g. coordinates in a PCA or MDS space).
1.  `mobest` provides a framework to perform spatiotemporal
    interpolation using Gaussian process regression (kriging) with the
    [`laGP`](https://CRAN.R-project.org/package=laGP) package to
    reconstruct an ancestry field based on the ancestry measure you
    provided.
2.  `mobest` allows to derive a similarity probability for samples of
    interest within the interpolated field, which – under certain
    circumstances – can be interpreted as an origin probability.
3.  `mobest` finally introduces functions to estimate and summarize a
    measure of mobility for the samples of interest, based on the
    similarity probability field.

Here is a simple, artificial example how 2. can be used:

``` r
library(magrittr)
set.seed(144)

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
  dependent = observations <- mobest::create_obs(
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
locate_product %>%
  ggplot() +
  geom_raster(
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  geom_point(
    mapping = aes(x = search_x, y = search_y), colour = "red"
  ) +
  coord_fixed() +
  ggtitle(paste0(
    "t for sample of interest = ", unique(locate_product$search_z), "\n",
    "t field time slice = ", unique(locate_product$field_z)
  ))
```

![](README_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Installation

Install the package from github with the following command in R:

    if(!require('remotes')) install.packages('remotes')
    remotes::install_github('nevrome/mobest')

## Overview

The following guide briefly lists the data types and functions in the
order you would usually call them and thus introduces the interface from
a technical point of view.

### Basic data types

`mobest` employs a number of basic [S3 data
types](http://adv-r.had.co.nz/S3.html) to formalize the input to almost
all of its functions. The constructors usually just check certain
properties to insure input correctness.

#### Spatial coordinates

`mobest::create_geopos` creates an object of class
`mobest_spatialpositions` which is a `tibble` that represents spatial
positions. Spatial positions in `mobest` are always 2-dimensional
coordinates in a Cartesian space. For real world coordinates that means,
that they have to be transformed to a projected coordinate system
(e.g. with `?sf::st_transform`): `mobest` can be used with longitude and
latitude coordinates.

``` r
mobest::create_geopos(
  id = 1:100,
  x = c(sample(100000:700000, 50), sample(300000:1000000, 50)),
  y = c(sample(100000:700000, 50), sample(300000:1000000, 50))
)
```

    ## # A tibble: 100 × 3
    ##       id      x      y
    ##    <int>  <int>  <int>
    ##  1     1 563893 248161
    ##  2     2 499058 241120
    ##  3     3 379661 597583
    ##  4     4 271319 396118
    ##  5     5 532892 614256
    ##  6     6 187987 398150
    ##  7     7 103322 395421
    ##  8     8 199717 464662
    ##  9     9 234251 102680
    ## 10    10 452401 352873
    ## # … with 90 more rows

For the interpolation fields we often want regular, spatial grids
covering a specific spatial area. These can be constructed with
`mobest::create_prediction_grid`, which takes an object of class `sf`
with polygons in a projected coordinate system. It also yields an object
of class `mobest_spatialpositions`.

Here is an example for the landmass of Europe, which we cover in a 250km
grid:

``` r
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

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

#### Spatiotemporal coordinates

`mobest_spatialpositions` can be transformed to
`mobest_spatiotemporalpositions` with `mobest::geopos_to_spatpos`. This
function calculates the permutations of all spatial positions with a
new, implicitly temporal dimension `z`:

``` r
mobest::geopos_to_spatpos(
  mobest::create_geopos(
    id = 1:100,
    x = c(sample(100000:700000, 50), sample(300000:1000000, 50)),
    y = c(sample(100000:700000, 50), sample(300000:1000000, 50))
  ),
  c(-5000, -4900, -4800)
)
```

    ## # A tibble: 300 × 5
    ##       id      x      y     z geo_id
    ##    <int>  <int>  <int> <dbl>  <int>
    ##  1     1 158229 205700 -5000      1
    ##  2     2 158229 205700 -4900      1
    ##  3     3 158229 205700 -4800      1
    ##  4     4 557976 121732 -5000      2
    ##  5     5 557976 121732 -4900      2
    ##  6     6 557976 121732 -4800      2
    ##  7     7 448428 527090 -5000      3
    ##  8     8 448428 527090 -4900      3
    ##  9     9 448428 527090 -4800      3
    ## 10    10 134001 354103 -5000      4
    ## # … with 290 more rows

`mobest::create_spatpos` directly creates an object of class
`mobest_spatiotemporalpositions` to represent spatiotemporal positions.

``` r
positions <- mobest::create_spatpos(
  id = 1:100,
  x = c(sample(100000:700000, 50), sample(300000:1000000, 50)),
  y = c(sample(100000:700000, 50), sample(300000:1000000, 50)),
  z = c(sample(-5000:-3500, 50), sample(-4500:-3000, 50))
)
```

    ## # A tibble: 100 × 4
    ##       id      x      y     z
    ##    <int>  <int>  <int> <int>
    ##  1     1 625599 599424 -3853
    ##  2     2 442475 155479 -3937
    ##  3     3 699770 539179 -4475
    ##  4     4 661104 254403 -4921
    ##  5     5 424049 640011 -3823
    ##  6     6 301860 528395 -4227
    ##  7     7 580030 109333 -4322
    ##  8     8 568890 269995 -3708
    ##  9     9 363356 128083 -3650
    ## 10    10 173665 250891 -4892
    ## # … with 90 more rows

#### Genetic coordinates

`mobest::create_obs` creates `mobest_observations`, which is a `tibble`
with genetic coordinates. Genetic coordinates can be any simple numeric
measure of ancestry, for example the position of the samples in PCA
space.

``` r
observations <- mobest::create_obs(
  ac1 = c(runif(50, 0, 0.6), runif(50, 0.4, 1)), # "ac" here for "ancestry component", e.g. PCA coordinate 1
  ac2 = c(runif(50, 0, 0.3), runif(50, 0.5, 1)) # e.g. PCA coordinate 2
)
```

    ## # A tibble: 100 × 2
    ##      ac1      ac2
    ##    <dbl>    <dbl>
    ##  1 0.107 0.219   
    ##  2 0.252 0.234   
    ##  3 0.368 0.000774
    ##  4 0.274 0.270   
    ##  5 0.188 0.300   
    ##  6 0.153 0.278   
    ##  7 0.480 0.0918  
    ##  8 0.113 0.152   
    ##  9 0.252 0.0133  
    ## 10 0.470 0.279   
    ## # … with 90 more rows

#### Kernel parameter settings

Gaussian process regression requires a parametrized covariance function:
a “kernel”. One `mobest_kernel` can be constructed with
`mobest::create_kernel`. `mobest_kernel` only represents one specific
kernel, though, for one specific genetic coordinate. Given that an
analysis typically involves multiple genetic coordinates (see
`mobest_observations`) `mobest::create_kernset` provides a wrapper to
bundle multiple kernels directly in an object of class
`mobest_kernelsetting`:

``` r
kernset <- mobest::create_kernset(
  ac1 = mobest::create_kernel(1000000, 1000000, 200, 0.1),
  ac2 = mobest::create_kernel(1000000, 1000000, 200, 0.1)
)
```

If a function requires both input of type `mobest_observations` and
`mobest_kernelsetting`, then the names of the individual ancestry
components must be identical.

#### Variability and permutations

When working with real data we often need to explore permutations of
data or account for uncertainty by sampling from distributions. To
represent that, `mobest` provides wrapper classes and constructors with
a `*_multi` suffix, to bundle multiple individual elements in a list
class. Some of the core functions provide interfaces that automatically
consider all permutations of these input lists.

Available are: - `mobest_spatialpositions_multi`
(`mobest::create_geopos_multi`) - `mobest_spatiotemporalpositions_multi`
(`mobest::create_spatpos_multi`) - `mobest_observations_multi`
(`mobest::create_obs_multi`) - `mobest_observationswitherror_multi`
(`mobest::create_obs_obserror_multi`) - `mobest_kernelsetting_multi`
(`mobest::create_kernset_multi`)

And here is an example how they can be filled with named arguments:

``` r
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

### Parameter estimation

One important question for the Gaussian process regression performed
within multiple of the core functions of `mobest` is a correct and
useful setting for the kernel parameters. The package therefore provides
different helper functions to estimate them.

#### Variogram calculation

`mobest::calculate_pairwise_distances` calculates different types of
pairwise distances (spatial, temporal, ancestry components) and returns
them in a long format `data.frame` object of class
`mobest_pairwisedistances`.

``` r
pairwise_distances <- mobest::calculate_pairwise_distances(
  independent = positions,
  dependent = observations,
  m_to_km = T
)
```

    ## # A tibble: 10,000 × 9
    ##      id1   id2 geo_dist time_dist obs_dist_total ac1_dist ac1_dist_resid
    ##    <int> <int>    <dbl>     <dbl>          <dbl>    <dbl>          <dbl>
    ##  1     1     1      0           0         0       0               0     
    ##  2     2     1    480.         84         0.146   0.146           0.332 
    ##  3     3     1     95.6       622         0.340   0.261           0.350 
    ##  4     4     1    347.       1068         0.175   0.167           0.415 
    ##  5     5     1    206.         30         0.114   0.0809          0.0858
    ##  6     6     1    331.        374         0.0755  0.0464          0.154 
    ##  7     7     1    492.        469         0.394   0.373           0.608 
    ##  8     8     1    334.        145         0.0668  0.00646         0.111 
    ##  9     9     1    539.        203         0.252   0.145           0.315 
    ## 10    10     1    571.       1039         0.368   0.363           0.665 
    ## # … with 9,990 more rows, and 2 more variables: ac2_dist <dbl>,
    ## #   ac2_dist_resid <dbl>

Helper functions are available to calculate the individual components of
this table:

``` r
geo_dist <- mobest::calculate_geo_pairwise_distances(positions)
time_dist <- mobest::calculate_time_pairwise_distances(positions)
obs_dist <- mobest::calculate_dependent_pairwise_distances(positions$id, observations)
```

`mobest::bin_pairwise_distances` bins the pairwise differences in an
object of class `mobest_pairwisedistances` and calculates an empirical
variogram (class `mobest_empiricalvariogram`) from them.

``` r
variogram <- mobest::bin_pairwise_distances(
  pairwise_distances,
  geo_bin = 0.1, time_bin = 100
)
```

    ## # A tibble: 4,805 × 8
    ##    geo_dist_cut time_dist_cut obs_dist_total ac1_dist ac1_dist_resid ac2_dist
    ##           <dbl>         <dbl>          <dbl>    <dbl>          <dbl>    <dbl>
    ##  1         0.05            50       0        0              0        0       
    ##  2         3.15           250       0.0695   0.0586         0.0472   0.0109  
    ##  3         4.35            50       0.00841  0.000226       0.000510 0.00818 
    ##  4         7.55           550       0.0760   0.0560         0.0365   0.0200  
    ##  5        12.2             50       0.0160   0.00323        0.00314  0.0127  
    ##  6        12.8            450       0.343    0.140          0.174    0.203   
    ##  7        13.4             50       0.00318  0.000859       0.00139  0.00232 
    ##  8        13.6            250       0.314    0.00560        0.00826  0.309   
    ##  9        13.8           1050       0.304    0.173          0.103    0.131   
    ## 10        14.4            550       0.000637 0.000467       0.000442 0.000170
    ## # … with 4,795 more rows, and 2 more variables: ac2_dist_resid <dbl>, n <int>

This variogram can for example be used to estimate the nugget parameter
of the GPR kernel settings, by filtering for pairwise “genetic”
distances with very small spatial and temporal distances.

#### Maximum likelihood estimation

`mobest::laGP_mle_anisotropic` wraps around `laGP::mleGPsep` to perform
marginal maximum likelihood inference for anisotropic (separable)
Gaussian lengthscale and nugget parameters.

``` r
mleGPsep_out <- mobest::laGP_mle_anisotropic(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2,
  verb = 0
)
```

    ## # A tibble: 4 × 9
    ##   mle_method dependent_var_id   dsx   dsy    dt      g optimizer_iterat… message
    ##   <chr>      <chr>            <dbl> <dbl> <dbl>  <dbl>             <int> <chr>  
    ## 1 mleGPsep   ac1              1482.  272. 1420. 0.0784                41 CONVER…
    ## 2 mleGPsep   ac1              1482.  272. 1420. 0.0784                41 CONVER…
    ## 3 mleGPsep   ac2               336.  311. 1881. 0.116                 17 CONVER…
    ## 4 mleGPsep   ac2               336.  311. 1881. 0.116                 17 CONVER…
    ## # … with 1 more variable: converged <int>

`mobest::laGP_jmle_anisotropic` does the same, but for joint maximum
likelihood inference.

``` r
jmleGPsep_out <- mobest::laGP_jmle_anisotropic(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2,
  verb = 0
)
```

    ## # A tibble: 4 × 9
    ##   mle_method dependent_var_id   dsx   dsy    dt      g optimizer_iterat… message
    ##   <chr>      <chr>            <dbl> <dbl> <dbl>  <dbl>             <int> <chr>  
    ## 1 jmleGPsep  ac1              1550. 1496. 1665. 0.0745                74 <NA>   
    ## 2 jmleGPsep  ac1              1550. 1496. 1665. 0.0745                74 <NA>   
    ## 3 jmleGPsep  ac2               336.  311. 1881. 0.116                130 <NA>   
    ## 4 jmleGPsep  ac2               336.  311. 1881. 0.116                130 <NA>   
    ## # … with 1 more variable: converged <int>

`mobest::laGP_mle_sequence_isotropic_fixed_g` implements a very specific
approach, where the mle is performed under the assumption of an
isotropic system, but with a series of scaling factors to explore the
space-time-relation. The nugget term g is fixed.

``` r
mle_sequence <- mobest::laGP_mle_sequence_isotropic_fixed_g(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2,
  g = 0.1,
  space_time_scaling_factor_sequence = seq(0.1, 2, 0.1),
  verb = 0
)
```

    ## # A tibble: 80 × 10
    ##    iteration dependent_var_id scaling_factor scaling_factor_fr… scaling_factor_…
    ##        <int> <chr>                     <dbl> <fractinl>         <fct>           
    ##  1         1 ac1                         0.1 0.1                1/10            
    ##  2         1 ac1                         0.2 0.2                1/5             
    ##  3         1 ac1                         0.3 0.3                3/10            
    ##  4         1 ac1                         0.4 0.4                2/5             
    ##  5         1 ac1                         0.5 0.5                1/2             
    ##  6         1 ac1                         0.6 0.6                3/5             
    ##  7         1 ac1                         0.7 0.7                7/10            
    ##  8         1 ac1                         0.8 0.8                4/5             
    ##  9         1 ac1                         0.9 0.9                9/10            
    ## 10         1 ac1                         1   1.0                1               
    ## # … with 70 more rows, and 5 more variables: d <dbl>, l <dbl>,
    ## #   optimizer_iterations <int>, ds <dbl>, dt <dbl>

#### Crossvalidation

`mobest::crossvalidate` allows to tackle the parameter estimation
challenge with simple cross-validation across a grid of kernel function
parameters. Internally it employs `mobest::create_model_grid` and
`mobest::run_model_grid` (see below).

``` r
kernels_to_test <- expand.grid(
  ds = seq(100,200, 50)*1000,
  dt = seq(100,200, 50)
) %>% purrr::pmap(function(...) {
    row <- list(...)
    mobest::create_kernset(
      ac1 = mobest::create_kernel(row$ds, row$ds, row$dt, 0.065),
      ac2 = mobest::create_kernel(row$ds, row$ds, row$dt, 0.08)
    )
  }) %>%
  magrittr::set_names(paste("kernel", 1:length(.), sep = "_")) %>%
  do.call(mobest::create_kernset_multi, .)

interpol_comparison <- mobest::crossvalidate(
  independent = positions,
  dependent = observations,
  kernel = kernels_to_test,
  iterations = 2,
  groups = 10,
  quiet = T
)
```

    ## # A tibble: 3,600 × 18
    ##    independent_table_id dependent_setting_id dependent_var_id kernel_setting_id
    ##    <fct>                <fct>                <chr>            <fct>            
    ##  1 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ##  2 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ##  3 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ##  4 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ##  5 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ##  6 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ##  7 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ##  8 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ##  9 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ## 10 ind_crossval_run_1   obs_crossval_run_1   ac1              kernel_1         
    ## # … with 3,590 more rows, and 14 more variables: pred_grid_id <fct>,
    ## #   mixing_iteration <int>, dsx <dbl>, dsy <dbl>, dt <dbl>, g <dbl>, id <int>,
    ## #   x <int>, y <int>, z <int>, mean <dbl>, sd <dbl>, measured <dbl>,
    ## #   difference <dbl>

### Spatiotemporal interpolation

The spatiotemporal interpolation workflow consists of the creation of a
list of models and then running each element in this list. The direct
interpolation function `mobest:::interpolate` has a minimal interface
and is therefore kept internal.

`mobest::create_model_grid` creates an object of class
`mobest_modelgrid` which holds all permutations of input elements. Each
row equals one complete model definition with all parameters and input
data fully defined.

``` r
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

    ## # A tibble: 32 × 9
    ##    independent_table_id dependent_setting_id dependent_var_id kernel_setting_id
    ##    <fct>                <fct>                <chr>            <fct>            
    ##  1 dating_1             obs1                 ac1              kernel_1         
    ##  2 dating_2             obs1                 ac1              kernel_1         
    ##  3 dating_1             obs2                 ac1              kernel_1         
    ##  4 dating_2             obs2                 ac1              kernel_1         
    ##  5 dating_1             obs1                 ac2              kernel_1         
    ##  6 dating_2             obs1                 ac2              kernel_1         
    ##  7 dating_1             obs2                 ac2              kernel_1         
    ##  8 dating_2             obs2                 ac2              kernel_1         
    ##  9 dating_1             obs1                 ac1              kernel_2         
    ## 10 dating_2             obs1                 ac1              kernel_2         
    ## # … with 22 more rows, and 5 more variables: pred_grid_id <fct>,
    ## #   independent_table <mbst_sp_>, dependent_var <named list>,
    ## #   kernel_setting <named list>, pred_grid <mbst_sp_>

`mobest::run_model_grid` runs each model and returns an unnested table
of interpolation results for each prediction grid point and each model
parameter setting.

``` r
interpol_grid <- mobest::run_model_grid(model_grid, quiet = T)
```

    ## # A tibble: 19,200 × 15
    ##    independent_table_id dependent_setting_id dependent_var_id kernel_setting_id
    ##    <fct>                <fct>                <chr>            <fct>            
    ##  1 dating_1             obs1                 ac1              kernel_1         
    ##  2 dating_1             obs1                 ac1              kernel_1         
    ##  3 dating_1             obs1                 ac1              kernel_1         
    ##  4 dating_1             obs1                 ac1              kernel_1         
    ##  5 dating_1             obs1                 ac1              kernel_1         
    ##  6 dating_1             obs1                 ac1              kernel_1         
    ##  7 dating_1             obs1                 ac1              kernel_1         
    ##  8 dating_1             obs1                 ac1              kernel_1         
    ##  9 dating_1             obs1                 ac1              kernel_1         
    ## 10 dating_1             obs1                 ac1              kernel_1         
    ## # … with 19,190 more rows, and 11 more variables: pred_grid_id <fct>,
    ## #   dsx <dbl>, dsy <dbl>, dt <dbl>, g <dbl>, id <int>, x <dbl>, y <dbl>,
    ## #   z <dbl>, mean <dbl>, sd <dbl>

### Origin search

`mobest::locate` uses the spatiotemporal interpolation to calculate a
similarity probability between input a set of “search” samples and
arbitrary spatiotemporal positions. It requires the necessary reference
sample input to perform the interpolation, which internally employs
`mobest::create_model_grid` and `mobest::run_model_grid`. The search
then yields a similarity probability value for each grid cell for each
search sample in an object of class `mobest_locateoverview`.

``` r
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

The spatiotemporal probability grids `locate` calculates are per
ancestry component (as put in via `dependent`/`search_dependent`). To
multiply the grids, `mobest` provides
`mobest::multiply_dependent_probabilities`, which yields an object of
class `mobest_locateproduct`.

``` r
mobest::multiply_dependent_probabilities(locate_simple)
```

    ## # A tibble: 800 × 14
    ##    search_id search_x search_y search_z independent_table_id dependent_setting_…
    ##        <int>    <int>    <int>    <int> <chr>                <chr>              
    ##  1         1   625599   599424    -3853 i                    d                  
    ##  2         1   625599   599424    -3853 i                    d                  
    ##  3         1   625599   599424    -3853 i                    d                  
    ##  4         1   625599   599424    -3853 i                    d                  
    ##  5         1   625599   599424    -3853 i                    d                  
    ##  6         1   625599   599424    -3853 i                    d                  
    ##  7         1   625599   599424    -3853 i                    d                  
    ##  8         1   625599   599424    -3853 i                    d                  
    ##  9         1   625599   599424    -3853 i                    d                  
    ## 10         1   625599   599424    -3853 i                    d                  
    ## # … with 790 more rows, and 8 more variables: field_z <dbl>,
    ## #   kernel_setting_id <fct>, pred_grid_id <fct>, field_id <int>, field_x <dbl>,
    ## #   field_y <dbl>, field_geo_id <int>, probability <dbl>

`mobest::locate` is actually just a special, simplified interface to
`mobest::locate_multi`, which adds another level of complexity, by
allowing multiple input values for `independent`, `dependent`, `kernel`,
`search_independent` and `search_dependent`. The result will consider
all permutations of these input settings (`independent` and
`search_independent` as well as `dependent` and `search_dependent` have
to be congruent, though).

``` r
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

`mobest::locate_multi` produces many probability grids for each sample.
Even after `mobest::multiply_dependent_probabilities` merges the
per-ancestry component iterations, that still leaves many parameter
permutations. `mobest::sum_probabilities_per_group` is a convenient
function to sum these to one, simplified probability grid of class
`mobest_locatesum`.

``` r
mobest::sum_probabilities_per_group(locate_product)
```

    ## # A tibble: 400 × 9
    ##    search_id field_id search_z search_x search_y field_x field_y field_z
    ##        <int>    <int>    <int>    <int>    <int>   <dbl>   <dbl>   <dbl>
    ##  1         1        1    -3853   625599   599424  100000  100000   -3853
    ##  2         1        2    -3853   625599   599424  200000  100000   -3853
    ##  3         1        3    -3853   625599   599424  300000  100000   -3853
    ##  4         1        4    -3853   625599   599424  400000  100000   -3853
    ##  5         1        5    -3853   625599   599424  500000  100000   -3853
    ##  6         1        6    -3853   625599   599424  600000  100000   -3853
    ##  7         1        7    -3853   625599   599424  700000  100000   -3853
    ##  8         1        8    -3853   625599   599424  800000  100000   -3853
    ##  9         1        9    -3853   625599   599424  900000  100000   -3853
    ## 10         1       10    -3853   625599   599424 1000000  100000   -3853
    ## # … with 390 more rows, and 1 more variable: probability <dbl>

`sum_probabilities_per_group` also allows to maintain the the
permutation groups, in case a full summary is not desired:

``` r
mobest::sum_probabilities_per_group(locate_product, dependent_setting_id, kernel_setting_id)
```

    ## # A tibble: 1,600 × 11
    ##    dependent_setting_id kernel_setting_id search_id field_id search_z search_x
    ##    <chr>                <fct>                 <int>    <int>    <int>    <int>
    ##  1 obs1                 kernel_1                  1        1    -3853   625599
    ##  2 obs1                 kernel_1                  1        2    -3853   625599
    ##  3 obs1                 kernel_1                  1        3    -3853   625599
    ##  4 obs1                 kernel_1                  1        4    -3853   625599
    ##  5 obs1                 kernel_1                  1        5    -3853   625599
    ##  6 obs1                 kernel_1                  1        6    -3853   625599
    ##  7 obs1                 kernel_1                  1        7    -3853   625599
    ##  8 obs1                 kernel_1                  1        8    -3853   625599
    ##  9 obs1                 kernel_1                  1        9    -3853   625599
    ## 10 obs1                 kernel_1                  1       10    -3853   625599
    ## # … with 1,590 more rows, and 5 more variables: search_y <int>, field_x <dbl>,
    ## #   field_y <dbl>, field_z <dbl>, probability <dbl>

### Origin vectors

To derive a simple, sample-wise measure of mobility,
`mobest::determine_origin_vectors` constructs what we call “origin
vectors” from objects of class `mobest_locateproduct`. Each vector
connects the spatial point where a respective individual was buried with
the point of highest genetic similarity in a respective search field.
The output is of class `mobest_originvectors` and documents distance and
direction of the “origin vector”, which can thus serve as a proxy of
mobility.

``` r
origin_vectors <- mobest::determine_origin_vectors(locate_product)
```

    ## # A tibble: 4 × 20
    ##   search_id search_x search_y search_z independent_tab… dependent_setti… field_z
    ##       <int>    <int>    <int>    <int> <chr>            <chr>              <dbl>
    ## 1         1   625599   599424    -3853 dating2          obs1               -3853
    ## 2         2   442475   155479    -3937 dating2          obs1               -3937
    ## 3         3   699770   539179    -4475 dating2          obs1               -4475
    ## 4         4   661104   254403    -4921 dating2          obs1               -4921
    ## # … with 13 more variables: kernel_setting_id <fct>, pred_grid_id <fct>,
    ## #   field_id <int>, field_x <dbl>, field_y <dbl>, field_geo_id <int>,
    ## #   probability <dbl>, ov_x <dbl>, ov_y <dbl>, ov_dist <dbl>, ov_dist_sd <dbl>,
    ## #   ov_angle_deg <dbl>, ov_angle_cut <chr>

Just as `mobest::sum_probabilities_per_group`, this summary can be split
to maintain the permutation groups introduced above.

``` r
mobest::determine_origin_vectors(locate_product, kernel_setting_id)
```

    ## # A tibble: 8 × 20
    ##   search_id search_x search_y search_z independent_tab… dependent_setti… field_z
    ##       <int>    <int>    <int>    <int> <chr>            <chr>              <dbl>
    ## 1         1   625599   599424    -3853 dating2          obs1               -3853
    ## 2         2   442475   155479    -3937 dating2          obs1               -3937
    ## 3         3   699770   539179    -4475 dating2          obs1               -4475
    ## 4         4   661104   254403    -4921 dating2          obs1               -4921
    ## 5         1   625599   599424    -3853 dating2          obs1               -3853
    ## 6         2   442475   155479    -3937 dating2          obs1               -3937
    ## 7         3   699770   539179    -4475 dating2          obs1               -4475
    ## 8         4   661104   254403    -4921 dating2          obs1               -4921
    ## # … with 13 more variables: kernel_setting_id <fct>, pred_grid_id <fct>,
    ## #   field_id <int>, field_x <dbl>, field_y <dbl>, field_geo_id <int>,
    ## #   probability <dbl>, ov_x <dbl>, ov_y <dbl>, ov_dist <dbl>, ov_dist_sd <dbl>,
    ## #   ov_angle_deg <dbl>, ov_angle_cut <chr>

In a very final step of the pipeline supported by `mobest` we can
summarise origin vectors through time.
`mobest::summarize_origin_vectors` allows for a moving window summary.
It also supports the deliberate grouping – note that also additional
variables can be introduced, e.g. a (spatial) region attribution of the
search samples.

``` r
origin_vectors$region_id <- c(
  "A", "B", "A", "C"
)

origin_summary <- mobest::summarize_origin_vectors(
  origin_vectors,
  region_id,
  window_start = -5000,
  window_stop = -3000,
  window_width = 100,
  window_step = 10
)
```

    ## # A tibble: 573 × 7
    ##    region_id     z undirected_mean_spatial_dist… directed_mean_s… mean_angle_deg
    ##    <chr>     <dbl>                         <dbl>            <dbl> <lgl>         
    ##  1 A         -4950                            NA               NA NA            
    ##  2 A         -4940                            NA               NA NA            
    ##  3 A         -4930                            NA               NA NA            
    ##  4 A         -4920                            NA               NA NA            
    ##  5 A         -4910                            NA               NA NA            
    ##  6 A         -4900                            NA               NA NA            
    ##  7 A         -4890                            NA               NA NA            
    ##  8 A         -4880                            NA               NA NA            
    ##  9 A         -4870                            NA               NA NA            
    ## 10 A         -4860                            NA               NA NA            
    ## # … with 563 more rows, and 2 more variables: se_spatial_distance <dbl>,
    ## #   sd_spatial_distance <dbl>

Empty time ranges can be identified with `mobest::find_no_data_windows`,
which is a minor, but useful helper function for plotting.

``` r
mobest::find_no_data_windows(
  origin_summary,
  region_id
)
```

    ## # A tibble: 6 × 3
    ##   region_id min_date_not_covered max_date_not_covered
    ##   <chr>                    <dbl>                <dbl>
    ## 1 A                        -4960                -4520
    ## 2 A                        -4430                -3900
    ## 3 A                        -3810                -3040
    ## 4 B                        -4960                -3980
    ## 5 B                        -3890                -3040
    ## 6 C                        -4880                -3040
