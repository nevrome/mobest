# How to cite mobest

If you use mobest for a scientific publication, please cite {cite:p}`Schmid2023`:

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

In the **materials section** you should

- document the specific software version you used
- explain the input data, specifically
  - how the input data was prepared and selected
  - how the spatial and temporal positions for each input sample was defined
  - how the position in dependent variable space was obtained (e.g. PCA on genotype data)
- list the kernel parameter settings and explain how they were obtained

Naturally **all code and data used to generate the mobest results should be shared** with the publication to ensure computational reproducibility.
