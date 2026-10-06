# Get the diet composition of a mizerReef model

Extends
[`mizer::getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.html)
so that the diet agrees with what mizerReef does when it projects the
model. Predators blocked by refuge (`blocked_pred = TRUE`) only
encounter the prey that are not hidden in refuge, `getVulnerable() * n`.
mizer's own method does not know about refuge and would compute their
diet from the whole prey population, overstating how much they eat of
every group that uses refuge. Because
[`mizer::plotDiet()`](https://sizespectrum.org/mizer/reference/plotDiet.html)
calls `getDiet()`, it shows the refuge-aware diet too.

## Usage

``` r
# S3 method for class 'mizerReef'
getDiet(
  object,
  proportion = TRUE,
  n = initialN(object),
  n_pp = initialNResource(object),
  n_other = initialNOther(object),
  ...,
  t = 0
)

# S3 method for class 'mizerReefSim'
getDiet(object, proportion = TRUE, time_range, drop = FALSE, ...)
```

## Arguments

- object:

  A `mizerReef` params object or a `mizerReefSim` object

- proportion:

  If `TRUE` (default) the function returns the diet as a proportion of
  the total consumption rate. If `FALSE` it returns the consumption rate
  in grams per year.

- n:

  A matrix of species abundances (species x size). Defaults to the
  initial abundances.

- n_pp:

  A vector of the resource abundance by size. Defaults to the initial
  resource abundance.

- n_other:

  A list of abundances for other dynamical components, such as algae and
  detritus. Defaults to the initial values.

- ...:

  Passed on to
  [`mizer::getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.html).

- t:

  For a params object, the time (a single year) at which the refuge and
  the feeding level are calculated. It matters when refuge degradation
  is switched on (see
  [`setDegradation()`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md))
  or another extension makes rates depend on time. Give it by name. For
  a `mizerReefSim`, each saved time is used; choose times with
  `time_range`.

- time_range:

  The range of times for which to return the diet, for a `mizerReefSim`.
  Defaults to all saved times. For a params object, use `t` instead.

- drop:

  If `TRUE`, dimensions of length 1 are removed from the array returned
  for a `mizerReefSim`.

## Value

For a params object, an array (predator species x predator size x prey),
as returned by
[`mizer::getDiet()`](https://sizespectrum.org/mizer/reference/getDiet.html).
For a `mizerReefSim`, an array with an additional first dimension for
time.

## Details

The diet of predators that are not blocked by refuge comes from the
whole prey population, as in mizer. For every consumer, the part of the
encounter that is eaten uses the model's feeding level at the given
abundances, algae and detritus biomasses and time.

mizer computes the diet from the feeding kernel itself rather than
through the model's encounter rate. So changes that other extensions
make to the encounter rate, such as a temperature effect, are not
included in the diet (sizespectrum/mizer#613).

## See also

[`getVulnerable()`](https://cmbeese.github.io/mizerReef/reference/getVulnerable.md)
