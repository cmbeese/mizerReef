# Total mortality rate in the reef ecosystem model

This function replaces the usual
[`mizerMort()`](https://sizespectrum.org/mizer/reference/mizerMort.html)
function and returns the sum of the usual mortality and size-based
external/ senescence mortality

## Usage

``` r
reefMort(params, n, n_pp, n_other, t, f_mort, pred_mort, ...)
```

## Arguments

- params:

  A MizerParams object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components of the ecosystem

- t:

  The time for which to do the calculation (Not used by standard mizer
  rate functions but useful for extensions with time-dependent
  parameters.)

- f_mort:

  A two dimensional array (species x size) with the fishing mortality

- pred_mort:

  A two dimensional array (species x size) with the predation mortality

- ...:

  Unused

## Value

A named two dimensional array (species x size) with the total mortality
rates.

## See also

Other mizer rate functions:
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md),
[`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md),
[`reefFeedingLevel()`](https://cmbeese.github.io/mizerReef/reference/reefFeedingLevel.md),
[`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md),
[`reefRates()`](https://cmbeese.github.io/mizerReef/reference/reefRates.md),
[`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
