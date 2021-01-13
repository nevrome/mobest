
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
    ##  1     1 288941 394406 -4577
    ##  2     2 234057 693448 -4245
    ##  3     3 224021 349221 -4995
    ##  4     4 326317 181117 -3849
    ##  5     5 465208 119391 -4122
    ##  6     6 293626 408280 -4333
    ##  7     7 669691 298818 -4952
    ##  8     8 689448 626573 -3738
    ##  9     9 597689 159969 -3637
    ## 10    10 502857 315369 -4808
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
    ##  1     1 288941 394406 394445
    ##  2     2 234057 693448 693419
    ##  3     3 224021 349221 349247
    ##  4     4 326317 181117 181057
    ##  5     5 465208 119391 119298
    ##  6     6 293626 408280 408364
    ##  7     7 669691 298818 298769
    ##  8     8 689448 626573 626539
    ##  9     9 597689 159969 160036
    ## 10    10 502857 315369 315270
    ## # … with 90 more rows
    ## 
    ## $run_b
    ## # A tibble: 100 x 4
    ##       id      x      y      z
    ##    <int>  <int>  <int>  <int>
    ##  1     1 288941 394406 394373
    ##  2     2 234057 693448 693414
    ##  3     3 224021 349221 349131
    ##  4     4 326317 181117 181090
    ##  5     5 465208 119391 119383
    ##  6     6 293626 408280 408312
    ##  7     7 669691 298818 298805
    ##  8     8 689448 626573 626632
    ##  9     9 597689 159969 159957
    ## 10    10 502857 315369 315467
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
    ## [1] 0.4317524 0.2741321 0.3127586 0.1454742 0.0455165 0.2347771
    ## 
    ## $ac2
    ## [1] 0.12160917 0.02178958 0.28612759 0.06515513 0.12168758 0.18770049

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
    ##  1     1     1      0           0         0       0               0     
    ##  2     2     1    304.        332         0.187   0.158           0.270 
    ##  3     3     1     79.1       418         0.203   0.119           0.0426
    ##  4     4     1    217.        728         0.292   0.286           0.293 
    ##  5     5     1    327.        455         0.386   0.386           0.397 
    ##  6     6     1     14.6       244         0.208   0.197           0.226 
    ##  7     7     1    393.        375         0.429   0.424           0.491 
    ##  8     8     1    463.        839         0.0917  0.0887          0.208 
    ##  9     9     1    388.        940         0.182   0.164           0.0488
    ## 10    10     1    228.        231         0.153   0.00644         0.0340
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

    ## # A tibble: 4,807 x 8
    ##    geo_dist_cut time_dist_cut obs_dist_total ac1_dist ac1_dist_resid ac2_dist
    ##           <dbl>         <dbl>          <dbl>    <dbl>          <dbl>    <dbl>
    ##  1         0.05            50         0       0              0       0       
    ##  2         4.75          1050         0.298   0.0325         0.0128  0.265   
    ##  3         5.75            50         0.0107  0.00534        0.00522 0.00533 
    ##  4         5.85          1150         0.351   0.278          0.201   0.0731  
    ##  5         7.85           550         0.0181  0.0172         0.0265  0.000905
    ##  6        10.8           1050         0.0463  0.0271         0.0550  0.0193  
    ##  7        11.6            150         0.567   0.419          0.404   0.148   
    ##  8        11.8             50         0.403   0.0217         0.0220  0.382   
    ##  9        14.6            250         0.0216  0.0194         0.0254  0.00218 
    ## 10        17.4            250         0.0458  0.0451         0.0537  0.000691
    ## # … with 4,797 more rows, and 2 more variables: ac2_dist_resid <dbl>, n <int>

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
    ## 1 mleGPsep   ac1                989. 1444. 1822. 0.0648    42 CONVERGENCE…     0
    ## 2 mleGPsep   ac1                989. 1444. 1822. 0.0648    42 CONVERGENCE…     0
    ## 3 mleGPsep   ac2               1042. 1211. 1377. 0.0803    32 CONVERGENCE…     0
    ## 4 mleGPsep   ac2               1042. 1211. 1377. 0.0803    32 CONVERGENCE…     0

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
    ## 1 jmleGPsep  ac1                 989. 1444. 1822. 0.0648    87 NA        0
    ## 2 jmleGPsep  ac1                 989. 1444. 1822. 0.0648    87 NA        0
    ## 3 jmleGPsep  ac2                1041. 1211. 1377. 0.0803    69 NA        0
    ## 4 jmleGPsep  ac2                1041. 1211. 1377. 0.0803    69 NA        0

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
    ##  1    94                1 100000   100     1 ac1_dist         0.0109 
    ##  2    94                1 100000   100     1 ac2_dist         0.119  
    ##  3    52                1 100000   100     1 ac1_dist         0.297  
    ##  4    52                1 100000   100     1 ac2_dist         0.0214 
    ##  5    37                1 100000   100     1 ac1_dist         0.512  
    ##  6    37                1 100000   100     1 ac2_dist         0.187  
    ##  7    27                1 100000   100     1 ac1_dist         0.256  
    ##  8    27                1 100000   100     1 ac2_dist         0.136  
    ##  9    73                1 100000   100     1 ac1_dist         0.00332
    ## 10    73                1 100000   100     1 ac2_dist         0.111  
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
