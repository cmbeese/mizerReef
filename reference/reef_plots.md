# Description of mizerReef plotting functions

In addition to the plotting functions offered by the mizer package,
MizerReef provides new and extended plotting functions to visualize both
input parameters and the results of dynamic simulations for reef
ecosystem models. Several functions adapt or build on mizer's plotting
tools for reef-specific features.

## Details

Available plotting functions:

**Input parameter plots**

|  |  |
|----|----|
| Plot | Description |
| [`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md) | Plots the proportion of individuals (by length) that are protected from predators for each species (the refuge profile) |
| [`plotDegradationScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegradationScale.md) | Plots a heatmap of the degradation scaling for refuge density in degradation simulations across bleaching years and size bins. |
| [`plotDegScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegScale.md) | Plots a faceted heatmap comparing all three built-in degradation trajectories (rubble, algae, recovery) side by side. |
| [`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md) | Plots vulnerability to predation by size and proportion for each species, either at steady state or for a chosen simulation time step. |
| [`plotRefugeDensity()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeDensity.md) | Plots the density of refuge by size at steady state and through time. |

**Result plots**

|  |  |
|----|----|
| Plot | Description |
| [`mizer::plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.html) | Plots the biomass of species and unstructured components through time (mizerReef extends this via a `getBiomass.mizerReefSim` method). |
| [`mizer::plotSpectra2()`](https://sizespectrum.org/mizer/reference/plotSpectra2.html) | Compare two size spectra (e.g., models or scenarios) in one plot. |
| [`mizer::plotSpectraRelative()`](https://sizespectrum.org/mizer/reference/plotSpectraRelative.html) | Plots relative difference between two spectra. |
| [`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md) | Plots the change (percent or relative proportion) between two spectra. |
| [`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md) | Relative contribution of each species group to total abundance, biomass, and productivity. |
| [`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md) | Plots productivity for each species through time or at steady state. |
| [`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md) | Productivity of two models or two size ranges in one plot. |
| [`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md) | Plots the relative or percent change in productivity between two models or scenarios. |
| [`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md) | Plots total abundance for each species at steady state. |
| [`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md) | Plots total biomass for each species within a given size range at steady state. |
| [`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md) | Total biomass of two models or two size ranges in one plot. |
| [`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md) | Relative change in total biomass between two models or scenarios. |

All plotting functions return a ggplot object, which can be further
customized using the ggplot2 grammar. Most functions accept either a
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html) or
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
object as input. Species colors and line types are controlled by the
`linecolour` and `linetype` slots in
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html),
and can be changed by the user.

For group or legend naming, you can add a column `group_names` to your
`species_params` data frame. This column should contain a nicely
formatted character string for each group or species, and will be used
in plot legends and facet labels for improved readability.

Most plots allow selection of a subset of species via the `species`
argument. The order of species in legends matches the species parameter
data frame.

## Usage

Call these functions with a
[MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html) or
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
object containing simulation results. See individual function
documentation for details and examples.

## References

For the original mizer plotting functions and further details, see:
https://sizespectrum.org/mizer/reference/plotting_functions.html

## Author

cmbeese
