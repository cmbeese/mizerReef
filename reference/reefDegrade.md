# Scales the refuge density by a given value at set times

Allows for the gradual degradation of habitat structure following an
acute disturbance by decreasing the availability of refuge over time.

## Usage

``` r
reefDegrade(params, n, n_pp, n_other, t, ...)
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

- ...:

  Unused

## Value

A numeric vector of scaled refuge densities for each size bin

## Details

The degradation is controlled by a scaling matrix `deg_scale` where:

- Column 1 represents the bleaching year (initial impact)

- Columns 2+ represent years 1, 2, 3... post-bleaching

At each timestep, the scaling factors are multiplied by the previous
timestep's refuge density.
[`setDegradation()`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md)
can also configure an independent, optional post-bleaching algae
growth/capacity boost (see
[`getAlgaeBoost()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeBoost.md));
that boost is unrelated to refuge density and is not computed by this
function.

## See also

Other mizer rate functions:
[`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md),
[`reefFeedingLevel()`](https://cmbeese.github.io/mizerReef/reference/reefFeedingLevel.md),
[`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md),
[`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md),
[`reefRates()`](https://cmbeese.github.io/mizerReef/reference/reefRates.md),
[`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
