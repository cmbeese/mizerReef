# Match observed growth rates

This does the same as
[`mizer::matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.html)
but in addition also rescales the consumption rates of algae and
detritus.

## Usage

``` r
matchReefGrowth(params, species = NULL, keep = c("egg", "biomass", "number"))
```

## Arguments

- params:

  A MizerParams object

- species:

  The species to be affected. Optional. By default all species for which
  growth information is available will be affected. A vector of species
  names, or a numeric vector with the species indices, or a logical
  vector indicating for each species whether it is to be affected (TRUE)
  or not.

- keep:

  A string determining which quantity is to be kept constant. The
  choices are "egg" which keeps the egg density constant, "biomass"
  which keeps the total biomass of the species constant and "number"
  which keeps the total number of individuals constant.

## Value

A modified MizerParams object with rescaled rates and rescaled species
parameters `gamma`,`h`, `ks` and `k`.
