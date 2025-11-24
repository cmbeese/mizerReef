# Plot the vulnerability to predation of species by weight

When called with a MizerParams object the initial vulnerability is
plotted. The complement of refuge.

## Usage

``` r
plotVulnerable(
  object,
  species = NULL,
  all.sizes = FALSE,
  return_data = FALSE,
  ...
)

plotlyVulnerable(
  object,
  species = NULL,
  all.sizes = FALSE,
  return_data = FALSE,
  ...
)
```

## Arguments

- object:

  An object of class MizerParams

- species:

  The species to be selected. Optional. By default all species are
  selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- all.sizes:

  If TRUE, then feeding level is plotted also for sizes outside a
  species' size range. Default FALSE.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default value is
  FALSE.

- ...:

  unused

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame the vulnerability at each size

## See also

plotting_functions,
[`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md),
[`plotRefuge()`](https://cmbeese.github.io/mizerReef/reference/plotRefuge.md)

Other plotting functions:
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotBiomass.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md),
[`plotRefuge()`](https://cmbeese.github.io/mizerReef/reference/plotRefuge.md),
[`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md),
[`plotSpectraRelative()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraRelative.md),
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md)
