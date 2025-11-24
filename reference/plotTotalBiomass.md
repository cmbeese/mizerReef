# Plot the total fishable biomass for each Species Group at steady state

This functions creates a barplot with the biomass of each Species Group
within a size range. Usually used in conjunction with
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)
to check for decoupling.

## Usage

``` r
plotTotalBiomass(
  object,
  species = NULL,
  min_fishing_l = NULL,
  max_fishing_l = NULL,
  return_data = FALSE,
  ...
)

plotlyTotalBiomass(object, ...)
```

## Arguments

- object:

  An object of class MizerParams

- species:

  The species to be selected. Optional. By default all species are
  selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- min_fishing_l:

  The minimum length (cm) for biomass estimates. Defaults to smallest
  size.

- max_fishing_l:

  The maximum length (cm) of for biomass estimates. Defaults to max
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

[`plotBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotBiomass.md),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md)

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
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
