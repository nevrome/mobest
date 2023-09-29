This is the code for the mobest documentation website. It uses the [Sphinx documentation generator](https://www.sphinx-doc.org/en/master).

The distinction of different software versions is handled with git tags (as implemented in [sphinx-multiversion](https://holzhaus.github.io/sphinx-multiversion/master/index.html)).

To build the website you need the sphinx CLI tool and a number of extensions for it. `sphinx.def` defines an [apptainer](https://apptainer.org) image with the complete setup.

The image can be build with
```
apptainer build docs/sphinx.sif docs/sphinx.def
```

and the website can then be rendered with
```
apptainer exec docs/sphinx.sif sphinx-multiversion docs docs/_build/html
```
