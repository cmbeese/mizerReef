# Scale reef abundances

Multiplies the abundances of all or of selected species by given factors
and then retunes the reproductive efficiency accordingly.

## Usage

``` r
scaleReefAbundance(params, factor)
```

## Arguments

- params:

  A mizer params object

- factor:

  The factor by which the abundance of each species is multiplied. This
  can be specified in two ways:

  - A named numeric vector where the name indicates the species and the
    value gives the factor for that species. Only the named species are
    affected.

  - A number that gives the factor for all foreground species.

## Value

An object of type
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)

## Details

Does not run the system to steady state. For that you should call
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
explicitly afterwards.
