# Remove some species from the model

This extends mizer's
[`removeSpecies()`](https://sizespectrum.org/mizer/reference/removeSpecies.html)
to also remove the relevant row from the detritus and algae consumption
arrays `rho`.

## Usage

``` r
# S3 method for class 'mizerReef'
removeSpecies(params, species, ...)
```

## Arguments

- params:

  A `mizerReef` object

- species:

  The species to be removed. A vector of species names, or a numeric
  vector of species indices, or a logical vector indicating for each
  species whether it is to be removed (TRUE) or not.

- ...:

  Unused

## Value

A `mizerReef` object with fewer species.
