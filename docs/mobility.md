# Diachronic mobility estimation with origin vectors

## Summarising origin vectors

In a very final step of the pipeline supported by `mobest`, we can combine origin vectors to meaningful summaries. Which summary statistics turn out to be useful strongly depends on the research questions guiding a particular project, so it is well likely that the following functions are not appropriate for a given use case.

`mobest::pack_origin_vectors` takes an object of class `mobest_originvectors` and merges iterations of origin vectors that might have emerged from permutations in `mobest::locate_multi` into a single, mean (!) origin vector for each search individual.

```r
mobest::pack_origin_vectors(origin_vectors, independent_table_id)
packed_origin_vectors <- mobest::pack_origin_vectors(origin_vectors)
```

An important limitation of `pack_origin_vectors`: The standard deviation (`ov_dist_sd`) and standard error of the mean (`ov_dist_se`) are recalculated from the origin vectors in the input `mobest_originvectors` object. The `ov_dist_sd` and `ov_dist_se` available there are not taken into account (as of now).

Beyond merging individual vector iterations, `mobest::summarize_origin_vectors` allows for a moving window summary across (!) origin vectors of individual samples. It accepts both objects of class `mobest_originvectors` and `mobest_originvectorspacked`. It also supports the deliberate grouping available for all functions following `multiply_dependent_probabilities` -- note that this explicitly includes additional variables that can be introduced even at this point in the pipeline, e.g. a (spatial) region attribution of the search samples.

```r
packed_origin_vectors$region_id <- c(
  "A", "B", "A", "C"
)

origin_summary <- mobest::summarize_origin_vectors(
  packed_origin_vectors,
  region_id,
  window_start = -5000,
  window_stop = -3000,
  window_width = 100,
  window_step = 10,
  quiet = T
)
```

Empty (i.e. not sufficiently informed from data) time ranges in this moving window summary can be identified with `mobest::find_no_data_windows`, which is a minor, but useful helper function e.g. for plotting.

```r
mobest::find_no_data_windows(
  origin_summary,
  region_id
)
```
