This R package provides types and functions for spatiotemporal
interpolation of human genetic ancestry components, similarity search
and the calculation of a derived measure for **mob**ility
**est**imation.

```{figure} img/example_movie.gif
---
align: right
width: 380px
---
```

- `mobest` assumes you have a set of genetic samples with spatial (two coordinates in a projected reference system) and temporal positions (years BC/AD) for which you calculated a derived, numeric measure of genetic ancestry (e.g. coordinates in a PCA or MDS space). 
- `mobest` provides a framework to perform spatiotemporal interpolation using Gaussian process regression (kriging) with the [`laGP`](https://CRAN.R-project.org/package=laGP) package to reconstruct an ancestry field based on the ancestry measure you provided.
- `mobest` allows to derive a similarity probability for samples of interest within the interpolated field, which – under certain circumstances – can be interpreted as an origin probability. See the example GIF on the right.
- `mobest` finally introduces functions to estimate and summarize a measure of mobility for the samples of interest, based on the similarity probability field.

The workflow in mobest version 1.0.0 was specifically developed to support this research compendium: <https://github.com/nevrome/mobest.analysis.2022>, which in turn underpins {cite:p}`Schmid2023`. If you use mobest, please cite this paper.

```none
@article{Schmid2023,
  doi = {10.1073/pnas.2218375120},
  url = {https://doi.org/10.1073/pnas.2218375120},
  year = {2023},
  month = feb,
  publisher = {Proceedings of the National Academy of Sciences},
  volume = {120},
  number = {9},
  author = {Clemens Schmid and Stephan Schiffels},
  title = {Estimating human mobility in Holocene Western Eurasia with large-scale ancient genomic data},
  journal = {Proceedings of the National Academy of Sciences}
}
```
