This is the code for the mobest documentation website. It uses the [Sphinx documentation generator](https://www.sphinx-doc.org/en/master).

To build the website you need the sphinx CLI tool and a number of extensions for it. `sphinx.def` defines an [apptainer](https://apptainer.org) image with the complete setup.

The image can be build with
```
apptainer build docs/sphinx.sif docs/sphinx.def
```

For development the latest state of the documentation can be rendered with
```
apptainer exec docs/sphinx.sif sphinx-build docs docs/_build/html
```

For deployment it should be build with [sphinx-multiversion](https://holzhaus.github.io/sphinx-multiversion/master/index.html), which will build separate documentation folders for each git tag (considers only committed changes!)
```
apptainer exec docs/sphinx.sif sphinx-multiversion docs docs/_build/html
```


