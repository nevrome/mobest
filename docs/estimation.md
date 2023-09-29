# Parameter estimation for optimal ancestry interpolation

One important question for the Gaussian process regression performed within multiple of the core functions of `mobest` is a correct and useful setting for the kernel parameters. The package therefore provides different helper functions to either estimate them or prepare data products that can be used to estimate them.

## Variogram calculation

`mobest::calculate_pairwise_distances` calculates different types of pairwise distances (spatial, temporal, ancestry components) for each input sample pair and returns them in a long format `tibble` of class `mobest_pairwisedistances`.

```r
pairwise_distances <- mobest::calculate_pairwise_distances(
  independent = positions,
  dependent = observations,
  m_to_km = T
)
```

Helper functions are available to calculate the individual components of this table:

```r
geo_dist <- mobest::calculate_geo_pairwise_distances(positions)
time_dist <- mobest::calculate_time_pairwise_distances(positions)
obs_dist <- mobest::calculate_dependent_pairwise_distances(positions$id, observations)
```

`mobest::bin_pairwise_distances` bins the pairwise distances in an `mobest_pairwisedistances` object and calculates an empirical variogram (class `mobest_empiricalvariogram`) from them.

```r
variogram <- mobest::bin_pairwise_distances(
  pairwise_distances,
  geo_bin = 100, time_bin = 100
)
```

This variogram can for example be used to estimate the nugget parameter of the GPR kernel settings, by filtering for pairwise "genetic" distances with very small spatial and temporal distances.

## Maximum likelihood estimation

`mobest::laGP_mle_anisotropic` wraps around `laGP::mleGPsep` to perform marginal maximum likelihood inference for anisotropic (separable) Gaussian lengthscale and nugget parameters.

```r
mleGPsep_out <- mobest::laGP_mle_anisotropic(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2,
  verb = 0
)
```

`mobest::laGP_jmle_anisotropic` does the same, but for joint maximum likelihood inference.

```r
jmleGPsep_out <- mobest::laGP_jmle_anisotropic(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2,
  verb = 0
)
```

`mobest::laGP_mle_sequence_isotropic_fixed_g` implements a very specific approach, where the mle is performed under the assumption of an isotropic system, but with a series of scaling factors to explore the space-time-relation. The nugget term g is fixed.

```r
mle_sequence <- mobest::laGP_mle_sequence_isotropic_fixed_g(
  independent = dplyr::mutate(positions, x = x/1000, y = y/1000),
  dependent = observations,
  iterations = 2,
  g = c(ac1 = 0.1, ac2 = 0.1),
  space_time_scaling_factor_sequence = seq(0.1, 2, 0.1),
  verb = 0
)
```

## Crossvalidation

`mobest::crossvalidate` allows to tackle the parameter estimation challenge with simple crossvalidation across a grid of kernel function parameters. Internally it employs `mobest::create_model_grid` and `mobest::run_model_grid` (see below). Crossvalidation is computationally expensive, but in our experience the best method for the kernel parameter estimation.

```r
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
