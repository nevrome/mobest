# Installing mobest

## Install mobest directly as an R package

mobest is an R package and can be **installed directly from GitHub** with the following code on the R console:

```r
if(!require('remotes')) install.packages('remotes')
remotes::install_github('nevrome/mobest')
```

You can also install specific/older versions of mobest with the following syntax: `nevrome/mobest[@ref|#pull|@*release]`. For example to install the publication release version you can run `remotes::install_github('nevrome/mobest@1.0.0')`.

For any of this to work a number of **system libraries** (mostly for processing geospatial data) have to be installed on your system, primarily for one particular dependency of mobest: the **`sf` R package**. The following table includes the libraries and the names of the relevant packages in the package management systems of various Linux distributions and MacOS.

| System library                                                        | deb package<br>(Ubuntu/Debian) | rpm package<br>(Fedora/CentOS) | pkgbuild package<br>(Arch) | brew package<br>(MacOS) |
|-----------------------------------------------------------------------|----------------------------------|---------------------------------|-------------------------|----------------------|
| [GDAL](https://gdal.org)                                              | libgdal-dev                      | gdal                            | gdal                    | gdal                 |
| [GEOS](https://libgeos.org/)                                          | libgeos-dev<br>libgeos++-dev        | geos-devel                      | geos                    | geos                 |
| [PROJ](https://proj.org)                                              | libproj-dev                      | proj-devel<br>sqlite-devel         | proj                    | proj                 |
| [UDUNITS-2](https://www.unidata.ucar.edu/software/udunits/)           | libudunits2-dev                  | udunits                         | udunits                 | udunits              |

The `sf` package maintainers provide a good explanation how to install these: <https://r-spatial.github.io/sf/#installing>

## Create an apptainer image to run mobest

If installing system libraries is not possible (for example because you don't have root access) or desirable, then mobest can also be run through a virtualization layer. One option is the high performance computing container system [**apptainer**](https://apptainer.org) (formerly ["singularity"](https://apptainer.org/news/community-announcement-20211130)). To do this you can follow these steps:

1. **Install apptainer** on a system where you have root access by following these instructions: <https://apptainer.org/docs/admin/main/installation.html>
2. Create a new text file to **define a container**, e.g. `apptainer_mobest.def`. This file will specify a virtual machine with all necessary dependencies. Here is a possible, minimal configuration:

```none
Bootstrap: docker
From: rocker/geospatial:4.3.1

%post
  R --slave -e 'remotes::install_github(repo = "nevrome/mobest")'
```

We build the mobest container image on top of a predefined docker image that already includes R, the tidyverse and all relevant geospatial dependencies as provided by the [Rocker Project](https://rocker-project.org/images/versioned/rstudio.html). We then simply install the latest version of mobest.

3. **Build an image** `apptainer_mobest.sif` from this container definition file `apptainer_mobest.def` with

```bash
apptainer build apptainer_mobest.sif apptainer_mobest.def
```
This requires sudo rights. It will take a couple of minutes, because the base image has to be downloaded and mobest has to be build within the container. The resulting container will probably require about 1.5GB of storage space.

You can test `apptainer_mobest.sif` by running a simple computation with the version of R that is installed within it, e.g. with `apptainer exec apptainer_mobest.sif Rscript -e 1+1`.

If this works and if you are not already working there anyway, then you can copy the image to the computing system where you actually want to run the mobest analysis.

4. **Run the R script** that specifies your mobest analysis through the apptainer container with a modified version of this code:

```bash
apptainer exec \
  --bind=/path/to/your/analysis/directory \
  path/to/your/apptainer_mobest.sif \
  Rscript path/to/your/mobestRscript.R \
  /
```

This should execute `mobestRscript.R` just as running it directly through `Rscript`.

The `--bind` option is required to allow the container access to the relevant section of the file system where the R script reads and writes input and output data.

Note that you can also submit this through a job scheduler, if you have an **HPC system** at your disposal for your mobest analyses. Here is a basic example for a [SGE](https://docs.oracle.com/cd/E19279-01/820-3257-12/n1ge.html) script `mobest.sh` that could be submitted via the `qsub` command to run `mobestRscript.R` on a computing cluster:

```bash
#!/bin/bash
#
#$ -S /bin/bash  # defines bash as the shell for execution
#$ -N mobest     # name of the command that will be listed in the queue
#$ -cwd          # change to the current directory
#$ -j y          # join error and standard output in one file
#$ -o ~/log      # standard output file or directory
#$ -q archgen.q  # queue
#$ -pe smp 5     # use X CPU cores
#$ -l h_vmem=10G # request XGb of memory
#$ -V            # load personal profile

date
apptainer exec \
  --bind=/path/to/your/analysis/directory \
  path/to/your/apptainer_mobest.sif \
  Rscript path/to/your/mobestRscript.R \
  /
date
exit 0
```

See {ref}`An HPC crossvalidation setup for large lengthscale parameter spaces <estimation:an HPC crossvalidation setup for large lengthscale parameter spaces>` for an application where running mobest on an HPC is advisable.
