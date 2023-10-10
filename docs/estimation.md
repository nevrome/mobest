# Parameter estimation for optimal ancestry interpolation

One important question for the Gaussian process regression performed within multiple of the core functions of `mobest` is a correct and useful setting for the kernel parameters (see {ref}`Kernel parameter settings <basic:kernel parameter settings>` in the basic workflow description). Supplementary Text 2 of {cite:p}`Schmid2023` discusses this in detail. `mobest` provides different helper functions to either estimate the parameters or prepare data products that can be used to estimate them. Here we explain a practical way to estimate the nugget and lengthscale values.

For this tutorial we will use the data introduced and prepared in {doc}`A basic similarity search workflow <basic>`, specifically a `samples_projected.csv` table prepared in {ref}`Reading the input samples <basic:reading the input samples>`.

You can download a script with the main workflow explained below including the required test data here:

- [simple_similarity_search.R](data/simple_similarity_search.R)
- [samples_projected.csv](data/samples_projected.csv)

## Using a subset of the variogram to estimate the nugget parameter

We start by loading the input data - individual ancient DNA samples with their spatiotemporal and genetic position.

```r
samples_projected <- readr::read_csv("docs/data/samples_projected.csv")
# you have to replace "data/docs/" with the path to your copy of the file
```

### Determining pairwise distances

`mobest::calculate_pairwise_distances` allows to calculate different types of pairwise distances (spatial, temporal, dependent variables/ancestry components) for each input sample pair and returns them in a long format `tibble` of class `mobest_pairwisedistances`.

```r
distances_all <- mobest::calculate_pairwise_distances(
  independent = mobest::create_spatpos(
    id = samples_projected$Sample_ID,
    x = samples_projected$x,
    y = samples_projected$y,
    z = samples_projected$Date_BC_AD_Median
  ),
  dependent = mobest::create_obs(
    C1 = samples_projected$MDS_C1,
    C2 = samples_projected$MDS_C2
  )
)
```

Helper functions are available to calculate the individual types of distances, if this is desired.

```r
geo_dist  <- mobest::calculate_geo_pairwise_distances(positions)
time_dist <- mobest::calculate_time_pairwise_distances(positions)
obs_dist  <- mobest::calculate_dependent_pairwise_distances(positions$id, observations)
```

`mobest_pairwisedistances` includes the following columns/variables.

|Column         |Description |
|:--------------|:-----------|
|id1            ||
|id2            ||
|geo_dist       ||
|time_dist      ||
|obs_dist_total ||
|C\*_dist       ||
|C\*_dist_resid ||

Note that `mobest::calculate_pairwise_distances()` also calculates the distances in dependent variables space on the residuals of a linear model informed from the spatiotemporal positions. This behaviour can be turned off by setting `with_resid = FALSE`.

This table allows us to easily visualize and analyse the pairwise distance properties of our dataset, for example with a set of correlation plots.

```{figure} img/estimation/....png
...
```

### Summarizing distances in an empirical variogram

`mobest::bin_pairwise_distances` bins the pairwise distances in an `mobest_pairwisedistances` object and calculates an empirical variogram (class `mobest_empiricalvariogram`) from them.

```r
variogram <- mobest::bin_pairwise_distances(
  distances_all,
  geo_bin = 100, time_bin = 100
)
```

`mobest_empiricalvariogram` includes these columns/variables.

|Column         |Description |
|:--------------|:-----------|
|geo_dist_cut   ||
|time_dist_cut  ||
|obs_dist_total ||
|C\*_dist       ||
|C\*_dist_resid ||
|n              ||

This variogram can be visualized in various ways, one of which is a simple raster of the sample counts.

```{figure} img/estimation/....png
...
```

### Estimating the nugget parameter

This variogram can for example be used to estimate the nugget parameter of the GPR kernel settings, by filtering for pairwise "genetic" distances with very small spatial and temporal distances.

```r
d_all_long <- distances_all %>% tidyr::pivot_longer(
  cols = c(C1_mds_u_dist_resid, C2_mds_u_dist_resid, C3_mds_u_dist_resid),
  names_to = "dist_type", values_to = "dist_val"
) %>%
  dplyr::mutate(
    dist_type = dplyr::recode(
      dist_type,
      C1_mds_u_dist_resid = "C1_mds_u",
      C2_mds_u_dist_resid = "C2_mds_u",
      C3_mds_u_dist_resid = "C3_mds_u"
    )
  )
```

```r
lower_left_variogram <- d_all_long %>%
  # filter for small temporal and spatial pairwise distances
  dplyr::filter(time_dist < 50 & geo_dist < 50) %>%
  dplyr::filter(id1 != id2) %>%
  dplyr::mutate(
    # rescaling of the dist val to a relative proportion
    dist_val_adjusted = dplyr::case_when(
      dist_type == "C1_mds_u" ~ 0.5*(dist_val^2/stats::var(janno_final$C1_mds_u)),
      dist_type == "C2_mds_u" ~ 0.5*(dist_val^2/stats::var(janno_final$C2_mds_u)),
      dist_type == "C3_mds_u" ~ 0.5*(dist_val^2/stats::var(janno_final$C3_mds_u)),
    )
  )
```

```r
estimated_nuggets <- lower_left_variogram %>%
  dplyr::group_by(dist_type) %>%
  dplyr::summarise(nugget = mean(dist_val_adjusted, na.rm = T)) %>%
  dplyr::mutate(
    dependent_var_id = gsub("_dist", "", dist_type)
  )
```

```{figure} img/estimation/....png
...
```

## Finding optimal lengthscale parameters with crossvalidation

`mobest::crossvalidate` allows to tackle the parameter estimation challenge with simple crossvalidation across a grid of kernel function parameters. Internally it employs `mobest::create_model_grid` and `mobest::run_model_grid` (see below). Crossvalidation is computationally expensive, but in our experience the best method for the kernel parameter estimation.

### Basic setup

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

### Analyzing the crossvalidation results

```{figure} img/estimation/....png
...
```

### HPC setup for large lengthscale parameter spaces



<!--
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
-->