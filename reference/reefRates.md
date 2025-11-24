# Get all rates needed to project a mizerReef model

Calls other rate functions in sequence and collects the results in a
list.

## Usage

``` r
reefRates(params, n, n_pp, n_other, t = 0, effort, rates_fns, ...)
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

- e from `mizerEReproAndGrowth()`

- e_repro from `mizerERepro()`

- e_growth from `mizerEGrowth()`

- pred_rate from `mizerPredRate()`

- pred_mort from
  [`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md)

- sen_mort from
  [`reefSenMort()`](https://cmbeese.github.io/mizerReef/reference/reefSenMort.md)

- f_mort from `mizerFMort()`

- mort from
  [`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md)

- rdi from `mizerRDI()`

- rdd from `BevertonHoltRDD()`

- resource_mort from `mizerResourceMort()`

However you can replace any of these rate functions by your own rate
function if you wish, see `setRateFunction()` for details.

## See also

Other mizer rate functions:
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md),
[`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md),
[`reefFeedingLevel()`](https://cmbeese.github.io/mizerReef/reference/reefFeedingLevel.md),
[`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md),
[`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md),
[`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
