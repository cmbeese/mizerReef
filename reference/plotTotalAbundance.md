# Plot the total abundance for each species in a size range

This functions creates a barplot with the abundance of each species
within a given size range.

## Usage

``` r
plotTotalAbundance(
  object,
  stack = FALSE,
  species = NULL,
  min_fishing_l = NULL,
  max_fishing_l = NULL,
  return_data = FALSE,
  ...
)

plotlyTotalAbundance(object, ...)
```

## Arguments

- object:

  An object of class
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
  or [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html).
  If a
  [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html)
  object is provided, the abundance at the last time step is used.

- stack:

  A boolean value that determines whether bars are separated by species.
  Defaults to FALSE. If true, returns a stacked barplot with the total
  biomass for each species instead of individual bars for each species.
  Useful for comparison between steady states.

- species:

  The species to be selected. Optional. By default all species are
  selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- min_fishing_l:

  The minimum length (cm) for abundance estimates. Defaults to smallest
  size.

- max_fishing_l:

  The maximum length (cm) of for abundance estimates. Defaults to max
  length.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default value is
  FALSE.

- ...:

  unused

## Value

A ggplot2 object

## See also

[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.html),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md)

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
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
