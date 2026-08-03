# Get total predation mortality rate needed to project mizer reef model

Calculates the total predation mortality rate \\\mu\_{p,i}(w_p)\\ (in
units of 1/year) on each prey species by prey size:

## Usage

``` r
reefPredMort(params, n, n_pp, n_other, t, pred_rate, vulnerable = NULL, ...)
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

- pred_rate:

  A two dimensional array (predator species x predator size) with the
  predation rate

- vulnerable:

  Array (species x size) with the proportion of individuals that are not
  protected from predation by refuge

- ...:

  Unused

## Value

A two dimensional array (prey species x prey size) with the predation
mortality

## Details

\$\$\mu\_{p.i}(w_p) = \sum_j {\tt pred\\rate}\_j(w_p)\\ V\_{ji}(w_p)\\
\theta\_{ji}.\$\$

You would not usually call this function directly but instead use
[`getPredMort()`](https://sizespectrum.org/mizer/reference/getPredMort.html),
which then calls this function.

## See also

Other mizer rate functions:
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md),
[`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md),
[`reefFeedingLevel()`](https://cmbeese.github.io/mizerReef/reference/reefFeedingLevel.md),
[`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md),
[`reefRates()`](https://cmbeese.github.io/mizerReef/reference/reefRates.md),
[`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
