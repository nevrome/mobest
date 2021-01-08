
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
genetic ancestry (e.g. coordinates in a Multidimensional scaling space).
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
  id = 1:5,
  x = c(100, 500, 300, 150, 400),
  y = c(100, 500, 300, 250, 400),
  z = c(-5000, -4000, -4500, -4250, -4700)
)
```

    ## # A tibble: 5 x 4
    ##      id     x     y     z
    ##   <int> <dbl> <dbl> <dbl>
    ## 1     1   100   100 -5000
    ## 2     2   500   500 -4000
    ## 3     3   300   300 -4500
    ## 4     4   150   250 -4250
    ## 5     5   400   400 -4700

`create_spatpos_multi` creates a list of
`mobest_spatiotemporalpositions` objects. It’s meant to represent
positional uncertainty, by providing multiple sets of spatiotemporal
coordinates for points with identical IDs.

``` r
uncertain_positions <- mobest::create_spatpos_multi(
  id = 1:5,
  x = rep(list(c(100, 500, 300, 150, 400)), 2),
  y = rep(list(c(100, 500, 300, 250, 400)), 2),
  z = list(c(-5000, -4000, -4500, -4250, -4700), c(-5050, -3950, 4450, -4200, -4600)),
  it = c("run_a", "run_b")
)
```

    ## $run_a
    ## # A tibble: 5 x 4
    ##      id     x     y     z
    ##   <int> <dbl> <dbl> <dbl>
    ## 1     1   100   100 -5000
    ## 2     2   500   500 -4000
    ## 3     3   300   300 -4500
    ## 4     4   150   250 -4250
    ## 5     5   400   400 -4700
    ## 
    ## $run_b
    ## # A tibble: 5 x 4
    ##      id     x     y     z
    ##   <int> <dbl> <dbl> <dbl>
    ## 1     1   100   100 -5050
    ## 2     2   500   500 -3950
    ## 3     3   300   300  4450
    ## 4     4   150   250 -4200
    ## 5     5   400   400 -4600

`create_obs` creates named lists of observations vectors (class
`mobest_observations`) corresponding to the spatiotemporal positions
defined above.

``` r
observations <- mobest::create_obs(
  ac1 = c(1, 2, 3, 2, 2), # "ac" for "ancestry component"
  ac2 = c(2, 1, 4, 3, 4)
)
```

    ## $ac1
    ## [1] 1 2 3 2 2
    ## 
    ## $ac2
    ## [1] 2 1 4 3 4
    ## 
    ## attr(,"class")
    ## [1] "mobest_observations" "list"

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

    ## # A tibble: 25 x 9
    ##     Var1  Var2 geo_dist time_dist obs_dist_total ac1_dist ac1_dist_resid
    ##    <int> <int>    <dbl>     <dbl>          <dbl>    <dbl>          <dbl>
    ##  1     1     1    0             0           0           0       0.      
    ##  2     2     1    0.566      1000           1.41        1       1.11e-16
    ##  3     3     1    0.283       500           2.83        2       1.50e+ 0
    ##  4     4     1    0.158       750           1.41        1       5.00e- 1
    ##  5     5     1    0.424       300           2.24        1       5.00e- 1
    ##  6     1     2    0.566      1000           1.41        1       1.11e-16
    ##  7     2     2    0             0           0           0       0.      
    ##  8     3     2    0.283       500           3.16        1       1.50e+ 0
    ##  9     4     2    0.430       250           2           0       5.00e- 1
    ## 10     5     2    0.141       700           3           0       5.00e- 1
    ## # … with 15 more rows, and 2 more variables: ac2_dist <dbl>,
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

    ## # A tibble: 8 x 8
    ##   geo_dist_cut time_dist_cut obs_dist_total ac1_dist ac1_dist_resid ac2_dist
    ##          <dbl>         <dbl>          <dbl>    <dbl>          <dbl>    <dbl>
    ## 1         0.05            50           0       0           0.           0   
    ## 2         0.15           150           0.5     0.5         5.00e- 1     0   
    ## 3         0.15           250           1.      0.5         5.00e- 1     0.5 
    ## 4         0.15           650           4.5     0           1.25e- 1     4.5 
    ## 5         0.15           750           1.      0.5         1.25e- 1     0.5 
    ## 6         0.25           450           3.17    0.833       7.50e- 1     2.33
    ## 7         0.45           250           2.25    0.25        1.25e- 1     2   
    ## 8         0.55           950           1.      0.5         6.16e-33     0.5 
    ## # … with 2 more variables: ac2_dist_resid <dbl>, n <int>

#### Maximum likelihood estimation

#### Crossvalidation

### Spatiotemporal interpolation

### Mobility estimation
