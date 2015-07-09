# seegMBG
### Streamlined Model-based Geostatistics Functions for the SEEG Research Group

This repository will contain a number of miscellaneous functions to streamline model-based geostatistical analyses as applied in several SEEG projects.

#### Installation

You can install and load this package directly from github using the devtools package as follows:

```r
library(devtools)
install_github('SEEG-Oxford/seegMBG')
library(seegMBG)
```

Note that the package depends on INLA, which isn't available from CRAN.
If you don't already have INLA installed, you'll have to install it before `seegMBG`.
At the time of writing, the best way to install the latest stable version of INLA is like this:

```r
install.packages('INLA', repos = 'http://www.math.ntnu.no/inla/R/stable')
```

though see the [INLA website](http://www.r-inla.org/download) for more up to date details
