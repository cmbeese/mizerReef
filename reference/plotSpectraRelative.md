# Plot the relative difference or percent change between two spectra

This plots a measure of the relative difference between the steady state
spectra of two mizer objects. The user can choose how this difference is
calculated. Let the spectra of the two objects be represented as
\\N_1(w)\\ and \\N_2(w)\\.

## Usage

``` r
plotSpectraRelative(
  object1,
  object2,
  species = NULL,
  power,
  diff_method = "percent_change",
  ...
)

plotlySpectraRelative(object1, object2, diff_method, ...)
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

- diff_method:

  The method to calculate the relative change between models. If
  `percent.change`, the percent change is calculated relative to the
  value from object 1 with formula 100\*(new-old)/old. If `rel.diff` the
  relative difference is returned given by (new - old)/(old + new).

- ...:

  Parameters passed to `plotSpectra()`

## Value

A ggplot2 object

## Details

If `diff_method` is given as `percent_change`, this function plots the
percent change, given by \$\$ 100\*(N_2(w) - N_1(w)) / (N_1(w)).\$\$

If `diff_method` is given as `rel_diff` the difference is calculated
relative to their average, so \$\$2 (N_2(w) - N_1(w)) / (N_2(w) +
N_1(w)).\$\$

The individual spectra are calculated by the `plotSpectra()` function
which is passed all additional arguments you supply. So you can for
example determine a size range over which to average the simulation
results via the `time_range` argument. See `plotSpectra()` for more
options.

Note that it does not matter whether the relative difference is
calculated for the number density or the biomass density or the biomass
density in log weight because the factors of \\w\\ by which the
densities differ cancels out in the relative difference.

## See also

Other plotting functions:
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotBiomass.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md),
[`plotRefuge()`](https://cmbeese.github.io/mizerReef/reference/plotRefuge.md),
[`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md),
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
