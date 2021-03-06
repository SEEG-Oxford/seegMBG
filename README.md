[![Build Status](https://travis-ci.org/SEEG-Oxford/seegMBG.svg)](https://travis-ci.org/SEEG-Oxford/seegMBG)
![CRAN](http://www.r-pkg.org/badges/version/seegMBG)

# seegMBG
## Streamlined Model-based Geostatistics Functions for the SEEG Research Group

This repository will contain a number of miscellaneous functions to streamline model-based geostatistical analyses as applied in several SEEG projects.

### Installation

##### Dependencies

`seegMBG` depends on the INLA R package which isn't available on CRAN.
If you don't already have INLA installed, you'll have to install it before `seegMBG`.
At the time of writing, the best way to install the latest stable version of INLA is like this:

```r
install.packages('INLA', repos = 'http://www.math.ntnu.no/inla/R/stable')
```

though see the [INLA website](http://www.r-inla.org/download) for more up-to-date details.

This package also depends on the SEEG package `seegSDM` which can be installed from [its own GitHub repository](https://github.com/SEEG-Oxford/seegSDM) with the `devtools` R package:

```r
devtools::install_github('SEEG-Oxford/seegSDM')
```

##### seegMBG

You can then install this package with `devtools`:

```r
devtools::install_github('SEEG-Oxford/seegMBG')
```

