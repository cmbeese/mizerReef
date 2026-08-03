# Get the biomass of species and unstructured components through time

Extends
[`mizer::getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.html)
to also include the algae and detritus biomasses alongside the species
biomasses, so that
[`mizer::plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.html)
(which calls `getBiomass()` internally) shows them without needing its
own override.

## Usage

``` r
# S3 method for class 'mizerReefSim'
getBiomass(object, ...)
```

## Arguments

- object:

  A `mizerReefSim` object

- ...:

  Passed on to
  [`mizer::getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.html)

## Value

An
[mizer::ArrayTimeBySpecies](https://sizespectrum.org/mizer/reference/ArrayTimeBySpecies.html)
object (time x species/component)
