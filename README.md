
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
  x = sample(100000:500000, 100),
  y = sample(500000:999999, 100),
  z = sample(-5000:-3000, 100)
)
```

    ## # A tibble: 100 x 4
    ##       id      x      y     z
    ##    <int>  <int>  <int> <int>
    ##  1     1 494501 643635 -4130
    ##  2     2 321245 853345 -3921
    ##  3     3 280405 820998 -4546
    ##  4     4 229675 679294 -3637
    ##  5     5 306416 866653 -4217
    ##  6     6 371600 587803 -3954
    ##  7     7 477114 679509 -4340
    ##  8     8 363513 625192 -3939
    ##  9     9 483140 951035 -4410
    ## 10    10 445778 869684 -4079
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
    ##  1     1 494501 643635 643612
    ##  2     2 321245 853345 853371
    ##  3     3 280405 820998 821008
    ##  4     4 229675 679294 679234
    ##  5     5 306416 866653 866653
    ##  6     6 371600 587803 587837
    ##  7     7 477114 679509 679560
    ##  8     8 363513 625192 625126
    ##  9     9 483140 951035 951039
    ## 10    10 445778 869684 869671
    ## # … with 90 more rows
    ## 
    ## $run_b
    ## # A tibble: 100 x 4
    ##       id      x      y      z
    ##    <int>  <int>  <int>  <int>
    ##  1     1 494501 643635 643731
    ##  2     2 321245 853345 853288
    ##  3     3 280405 820998 820987
    ##  4     4 229675 679294 679268
    ##  5     5 306416 866653 866588
    ##  6     6 371600 587803 587857
    ##  7     7 477114 679509 679508
    ##  8     8 363513 625192 625234
    ##  9     9 483140 951035 950997
    ## 10    10 445778 869684 869620
    ## # … with 90 more rows

`create_obs` creates named lists of observations vectors (class
`mobest_observations`) corresponding to the spatiotemporal positions
defined above.

``` r
observations <- mobest::create_obs(
  ac1 = runif(100, 0, 1), # "ac" for "ancestry component"
  ac2 = runif(100, 0, 1)
)
```

The first ten observations:

    ## $ac1
    ## [1] 0.6384085 0.7622658 0.9102939 0.9375692 0.4311193 0.4443036
    ## 
    ## $ac2
    ## [1] 0.07592261 0.78487463 0.25787917 0.03202743 0.40308005 0.18283680

### Parameter estimation

Gaussian process regression requires a parametrized covariance function.
`mobest` provides helper functions to estimate the relevant parameters
from the input data.

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
    ##  1     1     1      0           0          0        0             0     
    ##  2     2     1    272.        209          0.720    0.124         0.0276
    ##  3     3     1    278.        416          0.327    0.272         0.152 
    ##  4     4     1    267.        493          0.302    0.299         0.128 
    ##  5     5     1    292.         87          0.387    0.207         0.349 
    ##  6     6     1    135.        176          0.222    0.194         0.251 
    ##  7     7     1     39.9       210          0.591    0.539         0.542 
    ##  8     8     1    132.        191          0.718    0.297         0.368 
    ##  9     9     1    308.        280          0.785    0.296         0.232 
    ## 10    10     1    231.         51          0.669    0.578         0.662 
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

    ## # A tibble: 4,756 x 8
    ##    geo_dist_cut time_dist_cut obs_dist_total ac1_dist ac1_dist_resid ac2_dist
    ##           <dbl>         <dbl>          <dbl>    <dbl>          <dbl>    <dbl>
    ##  1         0.05            50        0        0.          0          0       
    ##  2         4.05          1350        0.712    3.53e-1     0.438      0.359   
    ##  3         6.65           450        0.120    2.40e-2     0.0180     0.0956  
    ##  4         6.75           950        0.0240   2.07e-3     0.00835    0.0219  
    ##  5         7.45           750        0.00177  9.30e-4     0.00000833 0.000842
    ##  6        10.6           1250        0.0314   8.81e-5     0.00258    0.0313  
    ##  7        10.8            750        0.319    7.89e-2     0.0596     0.240   
    ##  8        10.8            250        0.171    7.32e-2     0.0807     0.0981  
    ##  9        10.8            350        0.00220  1.55e-3     0.00251    0.000646
    ## 10        11.0             50        0.128    8.15e-2     0.0834     0.0467  
    ## # … with 4,746 more rows, and 2 more variables: ac2_dist_resid <dbl>, n <int>

#### Maximum likelihood estimation

`mobest::laGP_mle_anisotropic` wraps around `laGP::mleGPsep` to perform
marginal maximum likelihood inference for anisotropic (separable)
Gaussian lengthscale and nugget parameters.

``` r
mleGPsep_out <- mobest::laGP_mle_anisotropic(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2
)
```

    ## # A tibble: 4 x 9
    ##   mle_method ancestry_compone…    dx    dy    dt      g   its msg           conv
    ##   <chr>      <chr>             <dbl> <dbl> <dbl>  <dbl> <int> <chr>        <int>
    ## 1 mleGPsep   ac1                533. 1353. 1541. 0.0805    45 CONVERGENCE…     0
    ## 2 mleGPsep   ac1                533. 1353. 1541. 0.0805    45 CONVERGENCE…     0
    ## 3 mleGPsep   ac2               1079. 1166. 1575. 0.0745    31 CONVERGENCE…     0
    ## 4 mleGPsep   ac2               1079. 1166. 1575. 0.0745    31 CONVERGENCE…     0

`mobest::laGP_jmle_anisotropic` wraps around `laGP::mleGPsep` to perform
joint maximum likelihood inference for anisotropic (separable) Gaussian
lengthscale and nugget parameters.

``` r
jmleGPsep_out <- mobest::laGP_jmle_anisotropic(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2
)
```

    ## # A tibble: 4 x 9
    ##   mle_method ancestry_component    dx    dy    dt      g   its msg    conv
    ##   <chr>      <chr>              <dbl> <dbl> <dbl>  <dbl> <int> <lgl> <int>
    ## 1 jmleGPsep  ac1                 533. 1353. 1541. 0.0805   113 NA        0
    ## 2 jmleGPsep  ac1                 533. 1353. 1541. 0.0805   113 NA        0
    ## 3 jmleGPsep  ac2                1079. 1166. 1575. 0.0745    81 NA        0
    ## 4 jmleGPsep  ac2                1079. 1166. 1575. 0.0745    81 NA        0

#### Crossvalidation

### Spatiotemporal interpolation

### Mobility estimation
