# Simple permutations

In the example in {doc}`A basic similarity search workflow <basic>` we performed the similarity search with `mobest::locate()` for only one parameter permutation, keeping everything constant except the two dependent variables. But as already laid out above in {ref}`The mobest_locateoverview table <basic:the mobest_locateoverview table>`, mobest can automatically consider more parameter permutations, the most basic of which are directly available in `locate()`.

Please note that all parameter permutations will be multiplied with all other permutations, causing the number of runs to grow rapidly. If you, for example, submit five time slices and five search samples, the number of runs will be $5*5=25$ times bigger than for one time slice and sample. The permutation mechanism is explained in more detail in {doc}`Advanced features of the mobest package <advanced>`

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
```

<details>
<summary>Code for this figure.</summary>

```r
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
</details>

```{figure} img/basic/search_map_two_timeslices.png
Similarity search map plot for two different time slices.
```

### Summarizing multiple searches in one figure

Plotting individual similarity probability map plots does not scale well, if the number of searches grows. `mobest::determine_origin_vectors` offers a way to derive a simple, summary statistic, by generating what we call "origin vectors" from objects of class `mobest_locateproduct`. Each vector connects the spatial point where a sample was found with the point of **highest genetic similarity** in the interpolated search field and its permutations. The output is of class `mobest_originvectors` and documents distance and direction of the "origin vector".

The origin vector summary can also be applied for individual parameter permutations, which is especially relevant for more complex application involving `mobest::locate_multi` (see {ref}`Similarity search with permutations <advanced:similarity search with permutations>`), but already comes in handy for two different search times.

```r
origin_vectors <- mobest::determine_origin_vectors(search_product, search_time)
```

The resulting object of type `mobest_originvectors` features one row for each vector with the following variables. Note that many variables here reflect mean values in this case, as `determine_origin_vectors()` can summarize information across parameter permutations.

|Column               |Description |
|:--------------------|:-----------|
|independent_table_id |Identifier of the spatiotemporal position permutation|
|dependent_setting_id |Identifier of the dependent variable space position permutation|
|dependent_var_id     |Identifier of the dependent variable|
|kernel_setting_id    |Identifier of the kernel setting permutation|
|pred_grid_id         |Identifier of the spatiotemporal prediction grid|
|field_id             |Identifier of the spatiotemporal prediction point|
|field_x              |Spatial x axis coordinate of the prediction point|
|field_y              |Spatial y axis coordinate of the prediction point|
|field_z              |Temporal coordinate (age) of the prediction point|
|field_geo_id         |Identifier of the spatial prediction point|
|search_id            |Identifier of the search sample|
|search_x             |Spatial x axis coordinate of the search sample|
|search_y             |Spatial y axis coordinate of the search sample|
|search_z             |Temporal coordinate (age) of the search sample|
|search_time          |Search time as provided by the user in `locate()`'s `search_time` argument|
|probability          |Probability density calculated in `locate()`|
|ov_x                 |Length of the origin vector in x direction|
|ov_y                 |Length of the origin vector in y direction|
|ov_dist              |Length of the origin vector in space (Euclidean distance)|
|ov_dist_se           |Standard error of the mean of all vector lengths|
|ov_dist_sd           |Standard deviation of all vector lengths|
|ov_angle_deg         |Direction of the origin vector as an angle in degree (0-360°)|

One way of making use of this object is by highlighting the points of maximal similarity probability in the map plot.

<details>
<summary>Code for this figure.</summary>

```r
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
  geom_point(
    data = origin_vectors,
    mapping = aes(field_x, field_y),
    fill = "orange", shape = 24
  ) +
  geom_segment(
    data = origin_vectors,
    mapping = aes(
      x = search_x, y = search_y,
      xend = field_x, yend = field_y
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red"
  ) +
  geom_label(
    data = origin_vectors,
    mapping = aes(
      x = (field_x + search_x)/2, y = (field_y + search_y)/2,
      label = paste0(round(ov_dist/1000, -2), "km")
    ),
    fill = "white", colour = "red", size = 2
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
</details>

```{figure} img/basic/search_map_two_timeslices_with_ovs.png
Similarity search map plot for two different time slices including arrows to the point of maximum similarity probability.
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
origin_vectors <- mobest::determine_origin_vectors(search_product)
```

<details>
<summary>Code for this figure.</summary>

```r
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
</details>

```{figure} img/basic/search_map_two_samples.png
Similarity search map plot for two different search samples.
```

With `origin_vectors` we can summarize the results for both samples in a single figure.

<details>
<summary>Code for this figure.</summary>

```r
ggplot() +
  geom_sf(
    data = research_land_outline_3035,
    fill = "grey", color = NA
  ) +
  geom_sf(
    data = research_area_3035,
    fill = NA, colour = "red",
    linetype = "solid", linewidth = 1
  ) +
  geom_point(
    data = origin_vectors,
    mapping = aes(search_x, search_y),
    colour = "red"
  ) +
  geom_point(
    data = origin_vectors,
    mapping = aes(field_x, field_y),
    fill = "orange", shape = 24
  ) +
  geom_segment(
    data = origin_vectors,
    mapping = aes(
      x = search_x, y = search_y,
      xend = field_x, yend = field_y
    ),
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "red"
  ) +
  theme_bw() +
  theme(
    axis.title = element_blank()
  )
```
</details>

```{figure} img/basic/search_map_two_samples_in_one_plot.png
Prototype of a figure that combines the independent search results for two samples in one figure.
```