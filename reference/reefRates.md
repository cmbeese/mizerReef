# Get all rates needed to project a mizerReef model

Calls other rate functions in sequence and collects the results in a
list.

## Usage

``` r
reefRates(params, n, n_pp, n_other, t = 0, effort, rates_fns, ...)
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

- effort:

  The effort for each fishing gear

- rates_fns:

  Named list of the functions to call to calculate the rates. Note that
  this list holds the functions themselves, not their names.

- ...:

  Unused

## Value

List of rates.

## Details

By default this function returns a list with the following components:

- predation vulnerability from
  [`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)

- encounter from
  [`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md)

- feeding level from
  [`reefFeedingLevel()`](https://cmbeese.github.io/mizerReef/reference/reefFeedingLevel.md)

- e from
  [`mizerEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/mizerEReproAndGrowth.html)

- e_repro from
  [`mizerERepro()`](https://sizespectrum.org/mizer/reference/mizerERepro.html)

- e_growth from
  [`mizerEGrowth()`](https://sizespectrum.org/mizer/reference/mizerEGrowth.html)

- pred_rate from
  [`mizerPredRate()`](https://sizespectrum.org/mizer/reference/mizerPredRate.html)

- pred_mort from
  [`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md)

- sen_mort from
  [`reefSenMort()`](https://cmbeese.github.io/mizerReef/reference/reefSenMort.md)

- f_mort from
  [`mizerFMort()`](https://sizespectrum.org/mizer/reference/mizerFMort.html)

- mort from
  [`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md)

- rdi from
  [`mizerRDI()`](https://sizespectrum.org/mizer/reference/mizerRDI.html)

- rdd from
  [`BevertonHoltRDD()`](https://sizespectrum.org/mizer/reference/BevertonHoltRDD.html)

- resource_mort from
  [`mizerResourceMort()`](https://sizespectrum.org/mizer/reference/mizerResourceMort.html)

However you can replace any of these rate functions by overriding the
relevant `project*.mizerReef` S3 method (e.g.
`projectEncounter.mizerReef`), following the
[`NextMethod()`](https://rdrr.io/r/base/UseMethod.html)-based
composition pattern used throughout this file so that other extension
packages keep working.

## See also

Other mizer rate functions:
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md),
[`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md),
[`reefFeedingLevel()`](https://cmbeese.github.io/mizerReef/reference/reefFeedingLevel.md),
[`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md),
[`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md),
[`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
