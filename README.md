
<!-- Rmd -> md -->

# mobest

This R package provides a pipeline for spatiotemporal interpolation of
human genetic ancestry components and a derived measure for **mob**ility
**est**imation. The workflow in version 1.0 was specifically developed
to support this research compendium:
<https://github.com/nevrome/mobest.analysis.2022>. For broader
applications the code would probably have to be adjusted. This is not
designed as a general-purpose package, but rather a structured
collection of functions for said paper.

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

`mobest::create_spatpos` creates an object of class
`mobest_spatiotemporalpositions` which is a `data.frame` that represents
spatiotemporal positions.

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

`mobest::create_spatpos_multi` creates a list of
`mobest_spatiotemporalpositions` objects. It’s meant to represent
positional uncertainty, by providing multiple sets of spatiotemporal
coordinates for points with identical IDs.

``` r
uncertain_positions <- mobest::create_spatpos_multi(
  run_a = positions %>% dplyr::mutate(z = z + sample(-100:100, 100)),
  run_b = positions %>% dplyr::mutate(z = z + sample(-100:100, 100))
)
```

    ## $run_a
    ## # A tibble: 100 × 4
    ##       id      x      y     z
    ##    <int>  <int>  <int> <int>
    ##  1     1 288941 394406 -4538
    ##  2     2 234057 693448 -4274
    ##  3     3 224021 349221 -4969
    ##  4     4 326317 181117 -3909
    ##  5     5 465208 119391 -4215
    ##  6     6 293626 408280 -4249
    ##  7     7 669691 298818 -5001
    ##  8     8 689448 626573 -3772
    ##  9     9 597689 159969 -3570
    ## 10    10 502857 315369 -4907
    ## # … with 90 more rows
    ## 
    ## $run_b
    ## # A tibble: 100 × 4
    ##       id      x      y     z
    ##    <int>  <int>  <int> <int>
    ##  1     1 288941 394406 -4610
    ##  2     2 234057 693448 -4279
    ##  3     3 224021 349221 -5085
    ##  4     4 326317 181117 -3876
    ##  5     5 465208 119391 -4130
    ##  6     6 293626 408280 -4301
    ##  7     7 669691 298818 -4965
    ##  8     8 689448 626573 -3679
    ##  9     9 597689 159969 -3649
    ## 10    10 502857 315369 -4710
    ## # … with 90 more rows

`mobest::create_obs` creates named lists of observations vectors (class
`mobest_observations`) corresponding to the spatiotemporal positions
defined above.

``` r
observations <- mobest::create_obs(
  ac1 = c(runif(50, 0, 0.6), runif(50, 0.4, 1)), # "ac" here for "ancestry component"
  ac2 = c(runif(50, 0, 0.3), runif(50, 0.5, 1))
)
```

The first ten observations:

    ## # A tibble: 100 × 2
    ##        ac1    ac2
    ##      <dbl>  <dbl>
    ##  1 0.432   0.122 
    ##  2 0.274   0.0218
    ##  3 0.313   0.286 
    ##  4 0.145   0.0652
    ##  5 0.0455  0.122 
    ##  6 0.235   0.188 
    ##  7 0.00786 0.189 
    ##  8 0.520   0.0985
    ##  9 0.596   0.0417
    ## 10 0.425   0.274 
    ## # … with 90 more rows

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

    ## # A tibble: 10,000 × 9
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

`mobest::bin_pairwise_distances` bins the pairwise differences in an
object of class `mobest_pairwisedistances` and calculates an empirical
variogram (class `mobest_empiricalvariogram`) from them.

``` r
variogram <- mobest::bin_pairwise_distances(
  pairwise_distances,
  geo_bin = 0.1, time_bin = 100
)
```

    ## # A tibble: 4,807 × 8
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

    ## # A tibble: 4 × 9
    ##   mle_method ancestry_component    dx    dy    dt      g   its msg          conv
    ##   <chr>      <chr>              <dbl> <dbl> <dbl>  <dbl> <int> <chr>       <int>
    ## 1 mleGPsep   ac1                 989. 1444. 1822. 0.0648    42 CONVERGENC…     0
    ## 2 mleGPsep   ac1                 989. 1444. 1822. 0.0648    42 CONVERGENC…     0
    ## 3 mleGPsep   ac2                1042. 1211. 1377. 0.0803    32 CONVERGENC…     0
    ## 4 mleGPsep   ac2                1042. 1211. 1377. 0.0803    32 CONVERGENC…     0

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
    ##   mle_method ancestry_component    dx    dy    dt      g   its msg    conv
    ##   <chr>      <chr>              <dbl> <dbl> <dbl>  <dbl> <int> <lgl> <int>
    ## 1 jmleGPsep  ac1                 989. 1444. 1822. 0.0648    87 NA        0
    ## 2 jmleGPsep  ac1                 989. 1444. 1822. 0.0648    87 NA        0
    ## 3 jmleGPsep  ac2                1041. 1211. 1377. 0.0803    69 NA        0
    ## 4 jmleGPsep  ac2                1041. 1211. 1377. 0.0803    69 NA        0

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
  space_time_scaling_factor_sequence = c(seq(0.1, 0.9, 0.1), 1, seq(2, 10, 1)),
  verb = 0
)
```

    ## # A tibble: 76 × 10
    ##    iteration ancestry_component scaling_factor scaling_factor_… scaling_factor_…
    ##        <int> <chr>                       <dbl> <fractinl>       <fct>           
    ##  1         1 ac1                           0.1 0.1              1/10            
    ##  2         1 ac1                           0.2 0.2              1/5             
    ##  3         1 ac1                           0.3 0.3              3/10            
    ##  4         1 ac1                           0.4 0.4              2/5             
    ##  5         1 ac1                           0.5 0.5              1/2             
    ##  6         1 ac1                           0.6 0.6              3/5             
    ##  7         1 ac1                           0.7 0.7              7/10            
    ##  8         1 ac1                           0.8 0.8              4/5             
    ##  9         1 ac1                           0.9 0.9              9/10            
    ## 10         1 ac1                           1   1.0              1               
    ## # … with 66 more rows, and 5 more variables: d <dbl>, l <dbl>, its <int>,
    ## #   ds <dbl>, dt <dbl>

#### Crossvalidation

`mobest::crossvalidate` allows to tackle the parameter estimation
challenge with simple cross-validation across a grid of kernel function
parameters. Internally it employs `mobest::create_model_grid` and
`mobest::run_model_grid` (see below).

``` r
par2try <- expand.grid(
  ds = seq(100,200, 50)*1000,
  dt = seq(100,200, 50)
)

kernels <- par2try %>% purrr::pmap(function(...) {
  row <- list(...)
  mobest::create_kernset(
    ac1 = mobest::create_kernel(row$ds, row$ds, row$dt, 0.065),
    ac2 = mobest::create_kernel(row$ds, row$ds, row$dt, 0.08)
  )
}) %>% magrittr::set_names(paste("kernel", par2try$ds/1000, par2try$dt, sep = "_"))

interpol_comparison <- mobest::crossvalidate(
  independent = positions,
  dependent = observations,
  kernel = kernels,
  iterations = 2,
  groups = 10
)
```

    ## # A tibble: 3,600 × 16
    ##    independent_table_id dependent_var_id kernel_setting_id kernel_dsx kernel_dsy
    ##    <fct>                <chr>            <fct>                  <dbl>      <dbl>
    ##  1 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ##  2 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ##  3 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ##  4 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ##  5 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ##  6 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ##  7 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ##  8 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ##  9 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ## 10 ind_crossval_run_1   ac1              kernel_100_100        100000     100000
    ## # … with 3,590 more rows, and 11 more variables: kernel_dt <dbl>,
    ## #   kernel_g <dbl>, id <int>, x <int>, y <int>, z <int>, mean <dbl>, sd <dbl>,
    ## #   mixing_iteration <int>, measured <dbl>, difference <dbl>

### Spatiotemporal interpolation

The spatiotemporal interpolation workflow consists of the creation of a
list of models and then running each element in this list.

`mobest::create_model_grid` creates an object of class
`mobest_modelgrid` which holds all permutations of input elements. Each
row equals one complete model definition with all parameters and input
data fully defined.

``` r
library(magrittr)
model_grid <- mobest::create_model_grid(
  independent = uncertain_positions,
  dependent = observations,
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
  prediction_grid = list(
    pred_grid = expand.grid(
      x = seq(100000, 1000000, 100000),
      y = seq(100000, 1000000, 100000),
      z = seq(-5500, -3000, 500)
    ) %>% {mobest::create_spatpos(
      id = 1:nrow(.),
      x = .$x,
      y = .$y,
      z = .$z
    )}
  )
)
```

    ## # A tibble: 8 × 8
    ##   independent_table_id dependent_var_id kernel_setting_id pred_grid_id
    ##   <fct>                <chr>            <fct>             <fct>       
    ## 1 run_a                ac1              kernel_1          pred_grid   
    ## 2 run_b                ac1              kernel_1          pred_grid   
    ## 3 run_a                ac2              kernel_1          pred_grid   
    ## 4 run_b                ac2              kernel_1          pred_grid   
    ## 5 run_a                ac1              kernel_2          pred_grid   
    ## 6 run_b                ac1              kernel_2          pred_grid   
    ## 7 run_a                ac2              kernel_2          pred_grid   
    ## 8 run_b                ac2              kernel_2          pred_grid   
    ## # … with 4 more variables: independent_table <named list>,
    ## #   dependent_var <named list>, kernel_setting <named list>,
    ## #   pred_grid <named list>

The helper function `mobest::prediction_grid_for_spatiotemporal_area`
can be used to construct a regular, spatiotemporal grid for the
`prediction_grid` argument of `create_model_grid`. It uses
`sf::st_make_grid` to create the spatial grid for a specified input
region.

`mobest::run_model_grid` runs each model and returns an unnested table
of interpolation results for each prediction grid point and each model
parameter setting.

``` r
interpol_grid <- mobest::run_model_grid(model_grid, quiet = T)
```

    ## # A tibble: 4,800 × 14
    ##    independent_table_… dependent_var_id kernel_setting_… pred_grid_id kernel_dsx
    ##    <fct>               <chr>            <fct>            <fct>             <dbl>
    ##  1 run_a               ac1              kernel_1         pred_grid       1000000
    ##  2 run_a               ac1              kernel_1         pred_grid       1000000
    ##  3 run_a               ac1              kernel_1         pred_grid       1000000
    ##  4 run_a               ac1              kernel_1         pred_grid       1000000
    ##  5 run_a               ac1              kernel_1         pred_grid       1000000
    ##  6 run_a               ac1              kernel_1         pred_grid       1000000
    ##  7 run_a               ac1              kernel_1         pred_grid       1000000
    ##  8 run_a               ac1              kernel_1         pred_grid       1000000
    ##  9 run_a               ac1              kernel_1         pred_grid       1000000
    ## 10 run_a               ac1              kernel_1         pred_grid       1000000
    ## # … with 4,790 more rows, and 9 more variables: kernel_dsy <dbl>,
    ## #   kernel_dt <dbl>, kernel_g <dbl>, id <int>, x <dbl>, y <dbl>, z <dbl>,
    ## #   mean <dbl>, sd <dbl>

### Mobility estimation

`mobest::search_spatial_origin` takes a number of “search points” (with
`independent` and `dependent` variables) for which the spatial origin
search should be performed, as well as an `interpol_grid`, which defines
the field in which the origin search should be performed. It returns a
data.frame with one row for each search point and interpol_grid setting.
Each row contains the specifics of the genetically closest origin point.

``` r
origin_grid <- mobest::search_spatial_origin(
  independent = uncertain_positions,
  dependent = observations,
  interpol_grid = interpol_grid,
  rearview_distance = 300
)
```

    ## running field setting 1 with search points run_a

    ## running field setting 1 with search points run_b

    ## running field setting 2 with search points run_a

    ## running field setting 2 with search points run_b

    ## running field setting 3 with search points run_a

    ## running field setting 3 with search points run_b

    ## running field setting 4 with search points run_a

    ## running field setting 4 with search points run_b

    ## # A tibble: 400 × 20
    ##    search_id search_x search_y search_z search_ac1 search_ac2 origin_id origin_x
    ##        <int>    <int>    <int>    <int>      <dbl>      <dbl>     <int>    <dbl>
    ##  1         1   288941   394406    -4538    0.432       0.122        174   400000
    ##  2         2   234057   693448    -4274    0.274       0.0218       214   400000
    ##  3         3   224021   349221    -4969    0.313       0.286         62   200000
    ##  4         4   326317   181117    -3909    0.145       0.0652       303   300000
    ##  5         5   465208   119391    -4215    0.0455      0.122        207   700000
    ##  6         6   293626   408280    -4249    0.235       0.188        217   700000
    ##  7         7   669691   298818    -5001    0.00786     0.189         41   100000
    ##  8         8   689448   626573    -3772    0.520       0.0985       331   100000
    ##  9         9   597689   159969    -3570    0.596       0.0417       331   100000
    ## 10        10   502857   315369    -4907    0.425       0.274        192   200000
    ## # … with 390 more rows, and 12 more variables: origin_y <dbl>, origin_z <dbl>,
    ## #   origin_mean_ac1 <dbl>, origin_mean_ac2 <dbl>, origin_sd_ac1 <dbl>,
    ## #   origin_sd_ac2 <dbl>, search_points_id <chr>, field_id <int>,
    ## #   field_independent_table_id <fct>, field_kernel_setting_id <fct>,
    ## #   spatial_distance <dbl>, angle_deg <dbl>

From this output multiple data products can be derived (e.g. for
plotting) with `mobest::average_origin_searchid`,
`mobest::average_origin_moving_window` and `mobest::no_data_windows`.
