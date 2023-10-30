# mizerReef

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This is an extension package for the mizer package (<https://sizespectrum.org/mizer/>) that makes it easy to set up a mizer model for tropical coral reef ecosystems. It includes unstructured algae and detritus resources in addition to plankton. Predator-prey interactions are mediated by predation refuge to simulate benthic complexity.

This package was developed to support the creation and exploration of a generic model for tropical coral reefs. 

## Installation

You can install the mizerReef package from GitHub with

```{r}
remotes::install_github("cmbeese/reefOld")
```

If this gives an error saying “there is no package called remotes” then you also need to do

```{r}
install.packages("remotes")
```

before trying again to install mizerReef.

You may be prompted to update some of your existing packages. The two packages that you should always update are `mizer` and `mizerExperimental`, because the mizerReef package will always be designed to work with the most recent version of these.




