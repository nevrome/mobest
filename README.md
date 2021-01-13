
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
    ##  1     1 698352 664510 -4011
    ##  2     2 105800 611404 -4000
    ##  3     3 561038 120020 -4740
    ##  4     4 648682 413471 -3528
    ##  5     5 146553 626499 -3771
    ##  6     6 530289 397801 -3553
    ##  7     7 375846 690459 -3891
    ##  8     8 433355 665063 -4952
    ##  9     9 603837 629155 -4233
    ## 10    10 400885 597692 -3755
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
    ##  1     1 698352 664510 664433
    ##  2     2 105800 611404 611353
    ##  3     3 561038 120020 119964
    ##  4     4 648682 413471 413481
    ##  5     5 146553 626499 626523
    ##  6     6 530289 397801 397739
    ##  7     7 375846 690459 690545
    ##  8     8 433355 665063 665026
    ##  9     9 603837 629155 629110
    ## 10    10 400885 597692 597608
    ## # … with 90 more rows
    ## 
    ## $run_b
    ## # A tibble: 100 x 4
    ##       id      x      y      z
    ##    <int>  <int>  <int>  <int>
    ##  1     1 698352 664510 664477
    ##  2     2 105800 611404 611446
    ##  3     3 561038 120020 120090
    ##  4     4 648682 413471 413376
    ##  5     5 146553 626499 626558
    ##  6     6 530289 397801 397829
    ##  7     7 375846 690459 690386
    ##  8     8 433355 665063 665035
    ##  9     9 603837 629155 629153
    ## 10    10 400885 597692 597648
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
    ## [1] 0.19441118 0.19973072 0.31343616 0.31585774 0.02458499 0.39233757
    ## 
    ## $ac2
    ## [1] 0.02308397 0.19581604 0.25442538 0.10800100 0.03816481 0.26900888

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
    ##  1     1     1       0          0         0       0               0     
    ##  2     2     1     595.        11         0.173   0.00532         0.255 
    ##  3     3     1     562.       729         0.260   0.119           0.400 
    ##  4     4     1     256.       483         0.148   0.121           0.168 
    ##  5     5     1     553.       240         0.170   0.170           0.0381
    ##  6     6     1     315.       458         0.316   0.198           0.298 
    ##  7     7     1     324.       120         0.186   0.144           0.253 
    ##  8     8     1     265.       941         0.400   0.399           0.594 
    ##  9     9     1     101.       222         0.0613  0.0593          0.128 
    ## 10    10     1     305.       256         0.166   0.134           0.248 
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

    ## # A tibble: 4,812 x 8
    ##    geo_dist_cut time_dist_cut obs_dist_total ac1_dist ac1_dist_resid ac2_dist
    ##           <dbl>         <dbl>          <dbl>    <dbl>          <dbl>    <dbl>
    ##  1         0.05            50        0       0             0          0      
    ##  2         4.95           350        0.0185  0.00350       0.00132    0.0150 
    ##  3         7.25            50        0.00702 0.000398      0.000524   0.00662
    ##  4        10.2           1350        0.0303  0.000494      0.00427    0.0298 
    ##  5        11.6            750        0.248   0.0324        0.0159     0.215  
    ##  6        12.2             50        0.0112  0.000152      0.0000473  0.0110 
    ##  7        13.6            850        0.170   0.0795        0.0535     0.0909 
    ##  8        14.4           1050        0.262   0.187         0.132      0.0750 
    ##  9        17.4            250        0.107   0.0903        0.102      0.0162 
    ## 10        17.8            550        0.0137  0.000700      0.000280   0.0130 
    ## # … with 4,802 more rows, and 2 more variables: ac2_dist_resid <dbl>, n <int>

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
    ## 1 mleGPsep   ac1               1561. 1400. 1607. 0.0670    38 CONVERGENCE…     0
    ## 2 mleGPsep   ac1               1561. 1400. 1607. 0.0670    38 CONVERGENCE…     0
    ## 3 mleGPsep   ac2               1482. 1382. 1510. 0.0939    37 CONVERGENCE…     0
    ## 4 mleGPsep   ac2               1482. 1382. 1510. 0.0939    37 CONVERGENCE…     0

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
    ## 1 jmleGPsep  ac1                1561. 1400. 1607. 0.0670    89 NA        0
    ## 2 jmleGPsep  ac1                1561. 1400. 1607. 0.0670    89 NA        0
    ## 3 jmleGPsep  ac2                1482. 1382. 1510. 0.0939    55 NA        0
    ## 4 jmleGPsep  ac2                1482. 1382. 1510. 0.0939    55 NA        0

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
paramters. Internally it already employs `mobest::create_model_grid` and
`mobest::run_model_grid`.

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
    ##  1    76                1 100000   100     1 ac1_dist          0.198 
    ##  2    76                1 100000   100     1 ac2_dist          0.251 
    ##  3    63                1 100000   100     1 ac1_dist          0.118 
    ##  4    63                1 100000   100     1 ac2_dist          0.189 
    ##  5    24                1 100000   100     1 ac1_dist          0.0593
    ##  6    24                1 100000   100     1 ac2_dist          0.0432
    ##  7    88                1 100000   100     1 ac1_dist          0.339 
    ##  8    88                1 100000   100     1 ac2_dist          0.0807
    ##  9    22                1 100000   100     1 ac1_dist          0.257 
    ## 10    22                1 100000   100     1 ac2_dist          0.0944
    ## # … with 3,590 more rows

### Spatiotemporal interpolation

### Mobility estimation
