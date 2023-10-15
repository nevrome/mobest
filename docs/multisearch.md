# Simple permutations

In the example above we performed the similarity search with `mobest::locate()` for only one parameter permutation, keeping everything constant except the two dependent variables. But as already laid out above in {ref}`The mobest_locateoverview table <basic:the mobest_locateoverview table>`, mobest can automatically consider more parameter permutations, the most basic of which are directly available in `locate()`.

Please note that all parameter permutations will be multiplied with all other permutations, causing the number of runs to grow rapidly. If you, for example, submit five time slices and five search samples, the number of runs will be $5*5=25$ times bigger than for one time slice and sample.

## Multiple search time slices

As explained in {ref}`Search positions <basic:search positions>` the `search_time` argument can take an integer vector of relative or absolute ages. That means we can run the search not just for one, but for arbitrarily many time slices at with one call to `locate`.

Here is an example with two time slices.

```r
search_result <- mobest::locate(
  independent        = ind,
  dependent          = dep,
  kernel             = kernset,
  search_independent = search_ind,
  search_dependent   = search_dep,
  search_space_grid  = spatial_pred_grid,
  search_time        = c(-6500, -5700),
  search_time_mode   = "absolute"
)
search_product <- mobest::multiply_dependent_probabilities(search_result)

ggplot() +
  geom_raster(
    data = search_product,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  scale_fill_viridis_c() +
  geom_sf(
    data = research_area_3035,
    fill = NA, colour = "red",
    linetype = "solid", linewidth = 1
  ) +
  geom_point(
    data = search_samples,
    mapping = aes(x, y),
    colour = "red"
  ) +
  ggtitle(
    label = "<Stuttgart> ~5250BC",
    subtitle = "Early Neolithic (Linear Pottery Culture) - Lazaridis et al. 2014"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  ) +
  guides(
    fill = guide_colourbar(title = "Similarity\nsearch\nprobability")
  ) +
  facet_wrap(
    ~search_time,
    labeller = \(variable, value) {
      paste0("Search time: ", abs(value), "BC")
    }
  )
```

```{figure} img/basic/search_map_two_timeslices.png
Similarity search map plot for two different time slices.
```

## Multiple search samples

We can also select multiple search samples and prepare the input data for `mobest::locate()` 

```r
search_samples <- samples_projected %>%
  dplyr::filter(
    Sample_ID %in% c("Stuttgart_published.DG", "I5411")
  )

search_ind <- mobest::create_spatpos(
  id = search_samples$Sample_ID,
  x  = search_samples$x,
  y  = search_samples$y,
  z  = search_samples$Date_BC_AD_Median
)
search_dep <- mobest::create_obs(
  C1 = search_samples$MDS_C1,
  C2 = search_samples$MDS_C2
)
```

Here we set the `search_time_mode` to `"relative"` to get a different, meaningful search time for both samples.

```r
search_result <- mobest::locate(
  independent        = ind,
  dependent          = dep,
  kernel             = kernset,
  search_independent = search_ind,
  search_dependent   = search_dep,
  search_space_grid  = spatial_pred_grid,
  search_time        = -700,
  search_time_mode   = "relative"
)
search_product <- mobest::multiply_dependent_probabilities(search_result)

ggplot() +
  geom_raster(
    data = search_product,
    mapping = aes(x = field_x, y = field_y, fill = probability)
  ) +
  scale_fill_viridis_c() +
  geom_sf(
    data = research_area_3035,
    fill = NA, colour = "red",
    linetype = "solid", linewidth = 1
  ) +
  geom_point(
    data = search_samples %>% dplyr::rename(search_id = Sample_ID),
    mapping = aes(x, y),
    colour = "red"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  ) +
  guides(
    fill = guide_colourbar(title = "Similarity\nsearch\nprobability")
  ) +
  facet_wrap(
    ~search_id,
    ncol = 2,
    labeller = labeller(
      search_id = c(
        "Stuttgart_published.DG" = paste(
          "<Stuttgart> ~5250BC",
          "Early Neolithic (Linear Pottery culture) - Lazaridis et al. 2014",
          "Search time: ~5950BC",
          sep = "\n"
        ),
        "I5411" = paste(
          "<I5411> ~6650BC",
          "Mesolithic (Iron Gates) - Mathieson et al. 2018",
          "Search time: ~7350BC",
          sep = "\n"
        )
      )
    )
  )
```

```{figure} img/basic/search_map_two_samples.png
Similarity search map plot for two different search samples.
```

<!--

### Summarizing multiple search samples

Plotting individual similarity probability map plots does not scale well, when one wants to 

To derive a simple, sample-wise measure of mobility, `mobest::determine_origin_vectors` constructs what we call "origin vectors" from objects of class `mobest_locateproduct`. Each vector connects the spatial point where a sample was found (so for ancient samples that is usually where the respective individual was buried) with the point of **highest genetic similarity** in the interpolated search field and its permutations. The output is of class `mobest_originvectors` and documents distance and direction of the "origin vector". Under certain circumstances this vector can serve as a proxy for mobility.

```r
mobest::determine_origin_vectors(locate_product, quiet = T)
```

Just as in `mobest::fold_probabilities_per_group`, the origin vector search can be performed per permutation of the groups introduced above.

```r
origin_vectors <- mobest::determine_origin_vectors(locate_product, independent_table_id, quiet = T)
```

-->