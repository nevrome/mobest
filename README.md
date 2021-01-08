# mobest

This R package provides a pipeline for spatiotemporal interpolation of human genetic ancestry components and a derived measure for **mob**ility **est**imation. The workflow in version 1.0 was specifically developed for this research paper:

> ..., Estimating human mobility in Holocene Western Eurasia with bulk ancient genetic data, ...

## Installation

Install the package from github with the following command in R:

```
if(!require('remotes')) install.packages('remotes')
remotes::install_github('nevrome/mobest')
```

## Workflow

`mobest` assumes you have a set of ancient DNA samples with spatial (two coordinates in a projected reference system) and temporal positions (years BC/AD) for which you calculated a derived, numeric measure of genetic ancestry (e.g. coordinates in a Multidimensional scaling space). This package now provides a framework to perform spatiotemporal interpolation using Gaussian process regression (kriging) with the [`laGP`](https://CRAN.R-project.org/package=laGP) package to reconstruct an ancestry field based on the ancestry measure you provided. `mobest` then allows to estimate a point-wise measure of mobility based on a search for “ancestry origin” positions with similar genetic make-up in the respective past.

The research paper linked above explains the details and background. The following guide just lists the functions in the order you would usually call them and introduces the interface from a technical point of view.

### Parameter estimation

### Spatiotemporal interpolation

### Mobility estimation

