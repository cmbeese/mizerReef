# Remove some species from the model

This calls
[`mizer::removeSpecies()`](https://sizespectrum.org/mizer/reference/removeSpecies.html)
and in addition removes the relevant row from the detritus and algae
consumption arrays `rho`.

## Usage

``` r
removeSpecies(params, species)
```

## Arguments

- params:

  A MizerParams object

- species:

  The species to be removed. A vector of species names, or a numeric
  vector of species indices, or a logical vector indicating for each
  species whether it is to be removed (TRUE) or not.

## Value

A MizerParams object with fewer species.
