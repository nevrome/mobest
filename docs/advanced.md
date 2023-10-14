# Advanced features of the mobest package

`mobest::locate()` uses spatiotemporal interpolation to calculate spatial similarity probability maps between a set of search samples and an interpolated ancestry field at specific time slices. This basic and arguably most important use case of the mobest package is documented in {doc}`A basic similarity search workflow <basic>`. `locate()` hides a lot of the complexity of mobest, though, and some of that is documented in the following sections.

## Gaussian process regression on top of a linear model

An important detail of mobest is that it typically performs the spatiotemporal interpolation for an individual dependent variable in a two-step process. It starts with a simple linear model that predicts the dependent variable based on the three independent space and time variables, thus detrending the entire dataset. The more advanced, complex and computationally expensive Gaussian process regression model is then only applied to the residuals of the linear models.

The detrending step can be turned off for individual variables in `mobest::create_kernel`, by setting `on_residuals = FALSE` (see {ref}`Kernel parameter settings <types:kernel parameter settings>`).

## Spatiotemporal interpolation permutations in a model grid

Spatiotemporal interpolation is the first main operation mobest undertakes. The actual interpolation is performed in an internal function `mobest:::interpolate()`, which has a minimal interface and does not have to be called directly by the user. This is not least, because the main workflow in `mobest::locate()` (see {doc}`A basic similarity search workflow <basic>`) generally requires multiple interpolation runs. mobest therefore offers a slightly more complex, two step interpolation workflow, which consists of the creation of a list of models and then subsequently running each element in this list to construct different ancestry fields. 

`mobest::create_model_grid` creates an object of class `mobest_modelgrid` which holds all permutations of the field-defining input objects. Each row equals one complete model definition with all parameters and input data fully defined. Here is this function is called with the necessary input data constructors. Note that unlike for `locate()` we here require the `*_multi` constructors to express iterations of all settings (for more about this see {ref}`Similarity search with permutations <advanced:similarity search with permutations>` below).

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

As already pointed out in {ref}`The mobest_locateoverview table <basic:the mobest_locateoverview table>`, `mobest::locate()` is a special, simplified interface for `mobest::locate_multi()`. This more general function adds another level of complexity, by allowing multiple input values for `independent`, `dependent`, `kernel`, `search_independent` and `search_dependent` through a set of data types and constructor functions with the suffix `*_multi` (see {ref}`Permutation data types <types:permutation data types>`). The result will consider all permutations of these input settings (`independent` and `search_independent` as well as `dependent` and `search_dependent` have to be congruent, though).

```r
library(magrittr)
search_result <- mobest::locate_multi(
  independent        = mobest::create_spatpos_multi(...),
  dependent          = mobest::create_obs_multi(...),
  kernel             = mobest::create_kernset_multi(...),
  prediction_grid    = mobest::create_spatpos(...),
  search_independent = mobest::create_spatpos_multi(...),
  search_dependent   = mobest::create_obs_multi(...),
  search_space_grid  = mobest::create_spatpos(...),
  search_time        = 0,
  search_time_mode   = "relative"
)

search_product <- mobest::multiply_dependent_probabilities(search_result)
```

The result `locate_multi()` is also an object of type `mobest_locateoverview`, just as for `locate`. But `locate_multi()` actually makes use of its columns `independent_table_id`, `dependent_setting_id`, `dependent_var_id` and `kernel_setting_id`.

So `mobest::locate_multi` produces potentially many probability grids for each sample. Even after the per-dependent variable iterations are merged with `mobest::multiply_dependent_probabilities`, that could still leave many parameter permutations. `mobest::fold_probabilities_per_group` is a convenient function to combine these to a single, merged probability grid of class `mobest_locatefold`. The folding operation can be set in the argument `folding_operation`, where the default is a simple sum. Again, in the default setting, the output probabilities are normalized per permutation.

```r
mobest::fold_probabilities_per_group(search_product)
```

`fold_probabilities_per_group` also allows to maintain all or some permutation groups, in case a full summary is not desired:

```r
mobest::fold_probabilities_per_group(search_product, dependent_setting_id, kernel_setting_id)
```

### Temporal resampling as permutation application

The introduction of `mobest::locate_multi()` above has been abstract and devoid of any hint to the great power that comes with its permutation machinery. Here we want to show one application that is very relevant for archeogenetic and macro-archaeological applications of mobest: Temporal resampling.

Samples extracted from archaeological contexts usually have no precise age information that could be pin-pointed to exactly one year. Instead they can either be linked to a certain age range through relative chronology based on stratigraphy or typology, or they are dated with biological, physical or chemical methods of dating like for example [radiocarbon dating](https://en.wikipedia.org/wiki/Radiocarbon_dating). No matter how ingenious the method of dating might be, including sophisticated chronological modelling based on various lines of evidence, the outcome for the individual sample will always be a probability distribution over a set of potential years. And this set can be surprisingly large (up to several hundred years), with highly asymmetric probability distributions.

As a consequence, spatiotemporal interpolation only based on the median age, as demonstrated in {doc}`A basic similarity search workflow <basic>` is of questionable accuracy. The following section introduces a modified version of this simple workflow, now featuring temporal resampling based on archaeological age ranges and radiocarbon dates.

We start again by downloading and sub-setting data from {cite:p}`Schmid2023`. Here we already perform some additional steps in advance to reduce the complexity a bit. This includes the transformation of the spatial coordinates and the parsing and restructuring of the uncalibrated radiocarbon dates.

<details>
<summary>Code to prepare the input data table.</summary>

```r
# download .zip archives with tables from https://doi.org/10.17605/OSF.IO/6UWM5
utils::download.file(
  url = "https://osf.io/download/kej4s/",
  destfile = "docs/data/pnas_tables.zip"
)
# extract the relevant tables
utils::unzip(
  "docs/data/pnas_tables.zip",
  files = c("Dataset_S1.csv", "Dataset_S2.csv"),
  exdir = "docs/data/"
)
# read data files
samples_context_raw <- readr::read_csv("docs/data/Dataset_S1.csv")
samples_genetic_space_raw <- readr::read_csv("docs/data/Dataset_S2.csv")
# join them by sample name
samples_raw <- dplyr::left_join(
  samples_context_raw,
  samples_genetic_space_raw,
  by = "Sample_ID"
)
# create a useful subsets of this table
samples_selected <- samples_raw %>%
  dplyr::select(
    Sample_ID,
    Latitude, Longitude,
    Date_BC_AD_Start, Date_BC_AD_Median, Date_BC_AD_Stop, Date_C14,
    MDS_C1 = C1_mds_u, MDS_C2 = C2_mds_u
  )
# transform the spatial coordinates to EPSG:3035
samples_projected <- samples_selected %>%
  sf::st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326) %>%
  sf::st_transform(crs = 3035) %>% 
  dplyr::mutate(
    x = sf::st_coordinates(.)[,1],
    y = sf::st_coordinates(.)[,2],
    .after = "Sample_ID"
  ) %>%
  sf::st_drop_geometry()
# restructure radiocarbon information
samples_advanced <- samples_projected %>%
  tibble::add_column(
    .,
    purrr::map_dfr(
      .$Date_C14,
      function(y) {
        tibble::tibble(
          C14_ages = paste(stringr::str_match_all(y, ":(.*?)±")[[1]][,2], collapse = ";"),
          C14_sds = paste(stringr::str_match_all(y, "±(.*?)\\)")[[1]][,2], collapse = ";")
        )
      }
    ),
    .after = "Date_C14"
  )
readr::write_csv(samples_advanced, file = "docs/data/samples_advanced.csv")
```

</details>
<br>

You do not have to run this and can instead download the example table `samples_advanced.csv`. Additional to the simple `Date_BC_AD_Median` column this table has a different set of variables to express age.

| Variable | Type | Description |
|----------|------|-------------|
| Date_BC_AD_Start | int | The start of the age range for this sample in years BC/AD |
| Date_BC_AD_Stop  | int | The end of the age range for this sample in years BC/AD |
| Date_C14         | chr | A list of radiocarbon dates for each sample<br>in the form `(labCode-labNumber:ageBP±sd);...` |
| C14_ages         | chr | A string list with only the BP ages of the C14 dates |
| C14_sds          | chr | A string list with only the standard deviations of the C14 dates |

With `samples_advanced.csv` we can write code to draw age samples and then use them for the similarity search. The script and the data can be downloaded here:

- [temporal_resampling_similarity_search.R](data/temporal_resampling_similarity_search.R)
- [samples_advanced.csv](data/samples_advanced.csv)

#### Drawing ages for each sample

As expected we start this analysis by loading the relevant dependencies and data.

```r
library(magrittr)
library(ggplot2)

samples_advanced <- readr::read_csv("docs/data/samples_advanced.csv")
```

The first step to draw age samples is then to determine per-year probability densities for each sample. For the ones without radiocarbon dates this density follows a simple, uniform distribution. We assign each year a fraction, depending on the total number of years between `Date_BC_AD_Start` and `Date_BC_AD_Stop`. We can encode this in a simple helper function `contextual_date_uniform()`, which returns a tibble with the years in a column `age` and the densities in a column `sum_dens`.

```r
contextual_date_uniform <- function(startbcad, stopbcad) {
  tibble::tibble(
    age = startbcad:stopbcad,
    sum_dens = 1/(length(startbcad:stopbcad))
  )
}
```

For the samples with radiocarbon ages this is more involved, because we have to calibrate the radiocarbon ages and determine the normalized sum of their densities. For the calibration we use the Bchron R package {cite}`Haslett2008`, which implements a simple, fast [calibration algorithm](https://en.wikipedia.org/wiki/Radiocarbon_calibration). We then sum the per-year densities in case there are multiple uncalibrated dates for a given sample, normalize them and fill age range gaps in the `Bchron::BchronCalibrate()` output (some years fall below its density threshold) with 0 to get a continous sequence of years and densities. The output of `radiocarbon_date_sumcal()` has the same structure as `contextual_date_uniform()`: a tibble with the years in a column `age` and the densities in a column `sum_dens`.

```r
radiocarbon_date_sumcal <- function(ages, sds, cal_curve) {
  bol <- 1950 # c14 reference zero
  raw_calibration_output <- Bchron::BchronCalibrate(
    ages      = ages,
    ageSds    = sds,
    calCurves = rep(cal_curve, length(ages))
  )
  density_tables <- purrr::map(
    raw_calibration_output,
    function(y) {
      tibble::tibble(
        age = as.integer(-y$ageGrid + bol),
        densities = y$densities
      )
    }
  )
  sum_density_table <- dplyr::bind_rows(density_tables) %>%
    dplyr::group_by(age) %>%
    dplyr::summarise(
      sum_dens = sum(densities)/length(density_tables),
      .groups = "drop"
    ) %>%
    dplyr::right_join(
      ., data.frame(age = min(.$age):max(.$age)), by = "age"
    ) %>%
    dplyr::mutate(sum_dens = tidyr::replace_na(sum_dens, 0)) %>%
    dplyr::arrange(age)
  return(sum_density_table)
}
```

With these helper functions we can modify `samples_advanced` to include a new list-column `Date_BC_AD_Prob`, that features the density tibbles for each sample. Note that we set the calibration curve for all samples to `cal_curve = "intcal20"`. This is a sensible default for Western Eurasia, but not necessarily for other parts of the world.

```r
samples_with_age_densities <- samples_advanced %>%
  dplyr::mutate(
    Date_BC_AD_Prob = purrr::pmap(
      list(Date_BC_AD_Start, Date_BC_AD_Stop, C14_ages, C14_sds),
      function(context_start, context_stop, c14bps, c14sds) {
        if (!is.na(c14bps)) {
          radiocarbon_date_sumcal(
            ages = as.numeric(strsplit(c14bps, ";")[[1]]),
            sds = as.numeric(strsplit(c14sds, ";")[[1]]),
            cal_curve = "intcal20"
          )
        } else {
          contextual_date_uniform(
            startbcad = context_start,
            stopbcad = context_stop
          )
        }
      }
    )
  )
```

```{warning}
The calibration-based age probabilities we generate in this script are an improvement over using just median ages. But they still equate to a massive simplification of the per-sample age information. Each individual sample could potentially be informed by a dedicated chronological model to make the derived pear-year probabilities much more accurate and precise. But such models are typically not available for large meta-datasets like the one required for spatiotemporal interpolation on a continental scale.
```

In a last step we can define the number of age resampling runs we want to apply and draw this number of random samples from the age distributions for each sample. For the example here we chose two, but for a real world application a larger number (>50) is recommended.

```r
age_resampling_runs <- 2

samples_with_age_samples <- samples_with_age_densities %>%
  dplyr::mutate(
    Date_BC_AD_Samples = purrr::map(
      Date_BC_AD_Prob, function(x) {
        sample(
          x    = x$age,
          size = age_resampling_runs,
          prob = x$sum_dens, replace = T
        )
      }
    )
  )
```

`samples_with_age_samples` now includes another list column `Date_BC_AD_Samples` in which each cell features a vector of two individual ages. These are two possible ages for a given sample. Depending on the precision of the input age information and the shape of the radiocarbon calibration curve in the relevant age range, the individual age samples are often hundreds of years apart. This highlights the relevance of this age resampling exercise.

#### Applying the similarity search

```r
```