![GitHub R package
version](https://img.shields.io/github/r-package/v/nevrome/mobest)
[![R-CMD-check](https://github.com/nevrome/mobest/actions/workflows/check-release.yaml/badge.svg)](https://github.com/nevrome/mobest/actions/workflows/check-release.yaml)

# mobest

This R package provides types and functions for spatiotemporal interpolation of human genetic ancestry components, probabilistic similarity search and the calculation of a derived measure for **mob**ility **est**imation. The workflow in mobest version 1.0.0 was specifically developed to support [this](https://github.com/nevrome/mobest.analysis.2022) research compendium for [this](https://doi.org/10.1073/pnas.2218375120) paper.

**See the documentation at https://nevrome.de/mobest**

<img align="right" src="man/figures/example_movie.gif" width = 380>

- `mobest` assumes you have a set of genetic samples with spatial (two coordinates in a projected reference system) and temporal positions (years BC/AD) for which you calculated a derived, numeric measure of genetic ancestry (e.g. coordinates in a PCA or MDS space).
- `mobest` provides a framework to perform spatiotemporal interpolation using Gaussian process regression (kriging) with the [`laGP`](https://CRAN.R-project.org/package=laGP) package to reconstruct an ancestry field based on the ancestry measure you provided.
- `mobest` finally allows to derive a similarity probability for samples of interest within the interpolated field, which – under certain circumstances – can be interpreted as an origin probability.