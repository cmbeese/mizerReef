# Plot refuge density through time

Plots the refuge density at each size bin through time as colored lines,
showing how refuge degrades over the course of a simulation. Uses
viridis color scale to represent temporal succession.

## Usage

``` r
plotRefugeDensity(sim, return_data = FALSE, ...)

plotlyRefugeDensity(sim, ...)
```

## Arguments

- sim:

  A [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html)
  object

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default value is
  FALSE.

- ...:

  Unused

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame with refuge density at each size and time

## See also

[plotting_functions](https://sizespectrum.org/mizer/reference/plotting_functions.html),
[`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md),
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)

Other plotting functions:
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotDegScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegScale.md),
[`plotDegradationScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegradationScale.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md),
[`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md),
[`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md),
[`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md),
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
