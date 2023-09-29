# Installing mobest

## Directly as an R package

mobest is an R package and can be installed directly from GitHub on the R console with

```r
if(!require('remotes')) install.packages('remotes')
remotes::install_github('nevrome/mobest')
```

For this too work a number of system libraries (mostly for processing geospatial data) have to be installed on your system, primarily for the `sf` package. The following table includes the libraries and the names of the relevant packages in the package management systems of various Linux distributions and MacOS.

| System library                                                        | deb package<br>(Ubuntu/Debian) | rpm package<br>(Fedora/CentOS) | pkgbuild package<br>(Arch) | brew package<br>(MacOS) |
|-----------------------------------------------------------------------|----------------------------------|---------------------------------|-------------------------|----------------------|
| [GDAL](https://gdal.org)                                              | libgdal-dev                      | gdal                            | gdal                    | gdal                 |
| [GEOS](https://libgeos.org/)                                          | libgeos-dev<br>libgeos++-dev        | geos-devel                      | geos                    | geos                 |
| [PROJ](https://proj.org)                                              | libproj-dev                      | proj-devel<br>sqlite-devel         | proj                    | proj                 |
| [UDUNITS-2](https://www.unidata.ucar.edu/software/udunits/)           | libudunits2-dev                  | udunits                         | udunits                 | udunits              |

The `sf` package maintainers provide a good explanation how to install these: <https://r-spatial.github.io/sf/#installing>

## With apptainer
