# Reef feeding level

This function replaces the usual
[`mizerFeedingLevel()`](https://sizespectrum.org/mizer/reference/mizerFeedingLevel.html)
function and returns the a feeding level of 0 for piscivores.

## Usage

``` r
reefFeedingLevel(params, n, n_pp, n_other, t, encounter, ...)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
  object

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

- encounter:

  A two dimensional array (predator species x predator size) with the
  encounter rate.

- ...:

  Unused

## Value

A two dimensional array (predator species x predator size) with the
feeding level.

## See also

Other mizer rate functions:
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md),
[`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md),
[`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md),
[`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md),
[`reefRates()`](https://cmbeese.github.io/mizerReef/reference/reefRates.md),
[`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
