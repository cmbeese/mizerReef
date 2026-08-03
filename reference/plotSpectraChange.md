# Plot the change between two spectra

This plots the change between the steady state spectra of two mizer
objects. Let the spectra of the two objects be represented as \\N_1(w)\\
and \\N_2(w)\\. This function plots \$\$ \frac{N_2(w) -
N_1(w)}{N_1(w)}\$\$ expressed as a percentage (multiplied by 100) when
`use_percent = TRUE` (the default), or as the raw relative proportion
when `use_percent = FALSE`.

## Usage

``` r
plotSpectraChange(
  object1,
  object2,
  species = NULL,
  power,
  use_percent = TRUE,
  return_data = FALSE,
  ...
)

plotlySpectraChange(object1, object2, ...)
```

## Arguments

- object1:

  An object of class MizerSim or MizerParams

- object2:

  An object of class MizerSim or MizerParams

- species:

  The species to be selected. Optional. By default all species are
  selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each species
  whether it is to be selected (TRUE) or not.

- power:

  The abundance is plotted as the number density times the weight raised
  to this power. The default power = 1 gives the biomass density,
  whereas power = 2 gives the biomass density with respect to
  logarithmic size bins.

- use_percent:

  Logical. If TRUE (default), the change is expressed as a percentage
  (e.g. 50 for a 50% increase). If FALSE, the raw relative proportion is
  plotted instead (e.g. 0.5).

- return_data:

  Logical. If TRUE, returns the data frame underlying the plot instead
  of the plot itself. Default FALSE.

- ...:

  Parameters passed to `plotSpectra()`

## Value

A ggplot2 object, or a data frame if `return_data = TRUE`.

## Details

For the difference calculated relative to the average of the two
spectra, \\2 (N_2(w) - N_1(w)) / (N_2(w) + N_1(w))\\, use mizer's own
[`mizer::plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.html),
which already dispatches correctly for `mizerReef` objects.

The individual spectra are calculated by the
[`mizer::plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.html)
function which is passed all additional arguments you supply. So you can
for example determine a size range over which to average the simulation
results via the `time_range` argument. See
[`mizer::plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.html)
for more options.

## See also

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
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
