# Plot the vulnerability to predation of species by weight

Plot the vulnerability to predation of species by weight

## Usage

``` r
plotVulnerable(
  object,
  species = NULL,
  all.sizes = FALSE,
  time_step = NULL,
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

  An object of class
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
  or [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html).
  If a
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html)
  object is provided, the abundance at the last time step is used.

- species:

  The species to be selected. Optional. By default all species are
  selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- all.sizes:

  If TRUE, then vulnerability is plotted also for sizes outside a
  species' size range. Default FALSE.

- time_step:

  If `object` is a
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html)
  object, this optional parameter specifies which time step to use for
  calculating vulnerability. Default is the last time step.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default value is
  FALSE.

- ...:

  unused

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame the vulnerability at each size

## Details

When called with a
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
object the initial vulnerability is plotted. The complement of refuge.

## See also

[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.html),
[`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md),
[`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md)

Other plotting functions:
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotDegScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegScale.md),
[`plotDegradationScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegradationScale.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md),
[`plotRefugeDensity()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeDensity.md),
[`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md),
[`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md),
[`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md),
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md)
