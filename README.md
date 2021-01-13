
<!-- Rmd -> md -->

# mobest

This R package provides a pipeline for spatiotemporal interpolation of
human genetic ancestry components and a derived measure for **mob**ility
**est**imation. The workflow in version 1.0 was specifically developed
for this research paper:

> …, Estimating human mobility in Holocene Western Eurasia with bulk
> ancient genetic data, …

## Installation

Install the package from github with the following command in R:

    if(!require('remotes')) install.packages('remotes')
    remotes::install_github('nevrome/mobest')

## Workflow

`mobest` assumes you have a set of ancient DNA samples with spatial (two
coordinates in a projected reference system) and temporal positions
(years BC/AD) for which you calculated a derived, numeric measure of
genetic ancestry (e.g. coordinates in a multidimensional scaling space).
This package now provides a framework to perform spatiotemporal
interpolation using Gaussian process regression (kriging) with the
[`laGP`](https://CRAN.R-project.org/package=laGP) package to reconstruct
an ancestry field based on the ancestry measure you provided. `mobest`
then allows to estimate a point-wise measure of mobility based on a
search for “ancestry origin” positions with similar genetic make-up in
the respective past.

The research paper linked above explains the details and background. The
following guide just briefly lists the functions in the order you would
usually call them and introduces the interface from a technical point of
view.

### Basic classes for the input data

`create_spatpos` creates an object of class
`mobest_spatiotemporalpositions` which is a `data.frame` that represents
spatiotemporal positions.

``` r
positions <- mobest::create_spatpos(
  id = 1:100,
  x = c(sample(100000:700000, 50), sample(300000:999999, 50)),
  y = c(sample(100000:700000, 50), sample(300000:999999, 50)),
  z = c(sample(-5000:-3500, 50), sample(-4500:-3000, 50))
)
```

    ## # A tibble: 100 x 4
    ##       id      x      y     z
    ##    <int>  <int>  <int> <int>
    ##  1     1 370782 340270 -3713
    ##  2     2 123574 161629 -4438
    ##  3     3 202676 680325 -4344
    ##  4     4 552661 129982 -4986
    ##  5     5 238622 257656 -4246
    ##  6     6 595993 632414 -3675
    ##  7     7 148000 354861 -3806
    ##  8     8 101899 266849 -4572
    ##  9     9 110172 636817 -4021
    ## 10    10 524248 109962 -3687
    ## # … with 90 more rows

`create_spatpos_multi` creates a list of
`mobest_spatiotemporalpositions` objects. It’s meant to represent
positional uncertainty, by providing multiple sets of spatiotemporal
coordinates for points with identical IDs.

``` r
uncertain_positions <- mobest::create_spatpos_multi(
  id = 1:100,
  x = list(positions$x, positions$x),
  y = list(positions$y, positions$y),
  z = list(positions$y + sample(-100:100, 100), positions$y + sample(-100:100, 100)),
  it = c("run_a", "run_b")
)
```

    ## $run_a
    ## # A tibble: 100 x 4
    ##       id      x      y      z
    ##    <int>  <int>  <int>  <int>
    ##  1     1 370782 340270 340315
    ##  2     2 123574 161629 161610
    ##  3     3 202676 680325 680303
    ##  4     4 552661 129982 129957
    ##  5     5 238622 257656 257566
    ##  6     6 595993 632414 632451
    ##  7     7 148000 354861 354865
    ##  8     8 101899 266849 266914
    ##  9     9 110172 636817 636804
    ## 10    10 524248 109962 109918
    ## # … with 90 more rows
    ## 
    ## $run_b
    ## # A tibble: 100 x 4
    ##       id      x      y      z
    ##    <int>  <int>  <int>  <int>
    ##  1     1 370782 340270 340195
    ##  2     2 123574 161629 161618
    ##  3     3 202676 680325 680330
    ##  4     4 552661 129982 129899
    ##  5     5 238622 257656 257623
    ##  6     6 595993 632414 632421
    ##  7     7 148000 354861 354873
    ##  8     8 101899 266849 266922
    ##  9     9 110172 636817 636910
    ## 10    10 524248 109962 109968
    ## # … with 90 more rows

`create_obs` creates named lists of observations vectors (class
`mobest_observations`) corresponding to the spatiotemporal positions
defined above.

``` r
observations <- mobest::create_obs(
  ac1 = c(runif(50, 0, 0.6), runif(50, 0.4, 1)), # "ac" for "ancestry component"
  ac2 = c(runif(50, 0, 0.3), runif(50, 0.5, 1))
)
```

The first ten observations:

    ## $ac1
    ## [1] 0.4418675 0.3615023 0.2445929 0.3631534 0.5086189 0.3353639
    ## 
    ## $ac2
    ## [1] 0.13304016 0.04639533 0.19029199 0.20987723 0.14290661 0.11814348

### Parameter estimation

Gaussian process regression requires a parametrized covariance function
(a “kernel”). `mobest` provides helper functions to estimate the
relevant parameters from the input data.

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

    ## # A tibble: 10,000 x 9
    ##     Var1  Var2 geo_dist time_dist obs_dist_total ac1_dist ac1_dist_resid
    ##    <int> <int>    <dbl>     <dbl>          <dbl>    <dbl>          <dbl>
    ##  1     1     1       0          0         0        0              0     
    ##  2     2     1     305.       725         0.118    0.0804         0.161 
    ##  3     3     1     379.       631         0.205    0.197          0.191 
    ##  4     4     1     278.      1273         0.110    0.0787         0.0280
    ##  5     5     1     156.       533         0.0675   0.0668         0.204 
    ##  6     6     1     369.        38         0.108    0.107          0.320 
    ##  7     7     1     223.        93         0.0588   0.0358         0.0681
    ##  8     8     1     279.       859         0.125    0.0762         0.301 
    ##  9     9     1     395.       308         0.179    0.145          0.181 
    ## 10    10     1     277.        26         0.170    0.120          0.107 
    ## # … with 9,990 more rows, and 2 more variables: ac2_dist <dbl>,
    ## #   ac2_dist_resid <dbl>

`mobest::bin_pairwise_distances` bins the pairwaise differences in an
object of class `mobest_pairwisedistances` and calculates an empirical
variogram (class `mobest_empiricalvariogram`) from them.

``` r
variogram <- mobest::bin_pairwise_distances(
  pairwise_distances,
  geo_bin = 0.1, time_bin = 100
)
```

    ## # A tibble: 4,823 x 8
    ##    geo_dist_cut time_dist_cut obs_dist_total ac1_dist ac1_dist_resid ac2_dist
    ##           <dbl>         <dbl>          <dbl>    <dbl>          <dbl>    <dbl>
    ##  1         0.05            50        0       0            0           0.     
    ##  2         3.55          1150        0.0718  0.0708       0.0389      1.03e-3
    ##  3         5.55           650        0.143   0.00419      0.000808    1.38e-1
    ##  4         8.55            50        0.00224 0.000449     0.000488    1.79e-3
    ##  5         8.65           250        0.0462  0.0461       0.0393      8.04e-5
    ##  6         9.65           450        0.159   0.0742       0.0564      8.47e-2
    ##  7        10.0            850        0.00472 0.000463     0.00577     4.25e-3
    ##  8        10.4            350        0.0132  0.00162      0.000509    1.15e-2
    ##  9        13.2            150        0.00124 0.000372     0.00000915  8.68e-4
    ## 10        13.4           1050        0.0215  0.00496      0.0191      1.66e-2
    ## # … with 4,813 more rows, and 2 more variables: ac2_dist_resid <dbl>, n <int>

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

    ## # A tibble: 4 x 9
    ##   mle_method ancestry_compone…    dx    dy    dt      g   its msg           conv
    ##   <chr>      <chr>             <dbl> <dbl> <dbl>  <dbl> <int> <chr>        <int>
    ## 1 mleGPsep   ac1                857. 1438. 1380. 0.0637    35 CONVERGENCE…     0
    ## 2 mleGPsep   ac1                857. 1438. 1380. 0.0637    35 CONVERGENCE…     0
    ## 3 mleGPsep   ac2               1352. 1423.  926. 0.102     30 CONVERGENCE…     0
    ## 4 mleGPsep   ac2               1352. 1423.  926. 0.102     30 CONVERGENCE…     0

`mobest::laGP_jmle_anisotropic` wraps around `laGP::mleGPsep` to perform
joint maximum likelihood inference for anisotropic (separable) Gaussian
lengthscale and nugget parameters.

``` r
jmleGPsep_out <- mobest::laGP_jmle_anisotropic(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2,
  verb = 0
)
```

    ## # A tibble: 4 x 9
    ##   mle_method ancestry_component    dx    dy    dt      g   its msg    conv
    ##   <chr>      <chr>              <dbl> <dbl> <dbl>  <dbl> <int> <lgl> <int>
    ## 1 jmleGPsep  ac1                 857. 1438. 1380. 0.0637    92 NA        0
    ## 2 jmleGPsep  ac1                 857. 1438. 1380. 0.0637    92 NA        0
    ## 3 jmleGPsep  ac2                1353. 1423.  926. 0.102     89 NA        0
    ## 4 jmleGPsep  ac2                1353. 1423.  926. 0.102     89 NA        0

`mobest::laGP_mle_sequence_isotropic_fixed_g` implements a very specific
approach where the mle is performed under the assumption of an isotropic
system, but with a series of scaling factors for the
space-time-relation. The nugget term g is fixed.

``` r
mle_sequence <- mobest::laGP_mle_sequence_isotropic_fixed_g(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2,
  g = 0.1,
  space_time_scaling_factor_sequence = c(seq(0.1, 0.9, 0.1), 1, seq(2, 10, 1)),
  verb = 0
)
```

    ## # A tibble: 76 x 10
    ##    iteration ancestry_compon… scaling_factor scaling_factor_… scaling_factor_…
    ##        <int> <chr>                     <dbl> <fractinl>       <fct>           
    ##  1         1 ac1                         0.1 0.1              1/10            
    ##  2         1 ac1                         0.2 0.2              1/5             
    ##  3         1 ac1                         0.3 0.3              3/10            
    ##  4         1 ac1                         0.4 0.4              2/5             
    ##  5         1 ac1                         0.5 0.5              1/2             
    ##  6         1 ac1                         0.6 0.6              3/5             
    ##  7         1 ac1                         0.7 0.7              7/10            
    ##  8         1 ac1                         0.8 0.8              4/5             
    ##  9         1 ac1                         0.9 0.9              9/10            
    ## 10         1 ac1                         1   1.0              1               
    ## # … with 66 more rows, and 5 more variables: d <dbl>, l <dbl>, its <int>,
    ## #   ds <dbl>, dt <dbl>

#### Crossvalidation

`mobest::crossvalidate` allows to tackle the parameter estimation
challenge with simple cross-validation across a grid of kernel function
parameters. Internally it already employs `mobest::create_model_grid`
and `mobest::run_model_grid`.

``` r
interpol_comparison <- mobest::crossvalidate(
  independent = positions,
  dependent = observations,
  kernel = mobest::create_kernset_cross(
    ds = seq(100,200, 50)*1000,
    dt = seq(100,200, 50), 
    g = 0.1
  ),
  iterations = 2,
  groups = 10
)
```

    ## # A tibble: 3,600 x 7
    ##       id mixing_iteration     ds    dt     g dependent_var difference
    ##    <int>            <int>  <int> <int> <int> <chr>              <dbl>
    ##  1    51                1 100000   100     1 ac1_dist          0.144 
    ##  2    51                1 100000   100     1 ac2_dist          0.161 
    ##  3    39                1 100000   100     1 ac1_dist          0.306 
    ##  4    39                1 100000   100     1 ac2_dist          0.0495
    ##  5    96                1 100000   100     1 ac1_dist          0.173 
    ##  6    96                1 100000   100     1 ac2_dist          0.191 
    ##  7    65                1 100000   100     1 ac1_dist          0.0496
    ##  8    65                1 100000   100     1 ac2_dist          0.297 
    ##  9    46                1 100000   100     1 ac1_dist          0.0290
    ## 10    46                1 100000   100     1 ac2_dist          0.213 
    ## # … with 3,590 more rows

### Spatiotemporal interpolation

The spatiotemporal interpolation workflow consists of the creation of a
list of models and then running each element on this list.

`mobest::create_model_grid` creates an object of class
`mobest_modelgrid` which holds all permutations of input elements. Each
row equals one complete model definition with all parameters and input
data fully defined.

``` r
model_grid <- mobest::create_model_grid(
  independent = uncertain_positions,
  dependent = observations,
  kernel = mobest::create_kernset_multi(
    d = list(c(100000, 100000, 200)),
    g = 0.1,
    it = "kernel_100000_200_01"
  ),
  prediction_grid = mobest::create_spatpos_multi(
    id = c(1,2,3),
    x = list(sample(300000:999999, 3)),
    y = list(sample(300000:999999, 3)),
    z = list(sample(-4500:-3000, 3)),
    it = "pred_grid_1"
  )
)
```

    ## # A tibble: 4 x 8
    ##   independent_tab… dependent_var_id kernel_setting_… pred_grid_id
    ##   <fct>            <fct>            <fct>            <fct>       
    ## 1 run_a            ac1              kernel_100000_2… pred_grid_1 
    ## 2 run_b            ac1              kernel_100000_2… pred_grid_1 
    ## 3 run_a            ac2              kernel_100000_2… pred_grid_1 
    ## 4 run_b            ac2              kernel_100000_2… pred_grid_1 
    ## # … with 4 more variables: independent_table <named list>,
    ## #   dependent_var <mbst_bsr>, kernel_setting <named list>, pred_grid <named
    ## #   list>

The helper function `mobest::prediction_grid_for_spatiotemporal_area`
can be used to construct a regular, spatiotemporal grid for the
`prediction_grid` argument. It uses `sf::st_make_grid` under the hood to
create the spatial grid for an input region.

`mobest::run_model_grid`

### Mobility estimation
