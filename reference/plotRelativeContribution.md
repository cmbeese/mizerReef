# Plot the relative contribution of each species group to total abundance, total biomass, and total productivity

The group abundances, biomasses, productivities are calculated by the
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
and
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)
functions. These are passed all additional arguments you supply. See
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)
and
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)
for more details.

## Usage

``` r
plotRelativeContribution(
  object,
  min_size = NULL,
  min_fishing_l = NULL,
  include_inverts = FALSE,
  return_data = FALSE,
  ...
)

plotlyRelativeContribution(object, ...)
```

## Arguments

- object:

  An object of class
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)

- min_size:

  parameters be passed to
  [`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md)
  and
  [`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md).
  The minimum length (cm) of individuals for biomass estimates. Defaults
  to smallest size in the model.

- min_fishing_l:

  parameters be passed to
  [`getProductivity()`](https://cmbeese.github.io/mizerReef/reference/getProductivity.md).
  The minimum length (cm) of fished individuals for productivity
  estimates. Defaults to 7 cm.

- include_inverts:

  A boolean value that determines whether the "inverts" species group is
  included. Default is FALSE, since invertebrate productivity is
  typically not relevant to fishing yield.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default value is
  FALSE.

- ...:

  Arguments passed on to
  [`plotTotalBiomass`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
  [`plotTotalAbundance`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
  [`plotProductivity`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)

  `species`

  :   The species to be selected. Optional. By default all species are
      selected. A vector of species names, or a numeric vector with the
      species indices, or a logical vector indicating for each species
      whether it is to be selected (TRUE) or not.

  `max_fishing_l`

  :   The maximum length (cm) of for biomass estimates. Defaults to max
      length.

  `stack`

  :   A boolean value that determines whether bars are separated by
      species. Defaults to FALSE. If true, returns a stacked barplot
      with the total biomass for each species instead of individual bars
      for each species. Useful for comparison between steady states.

  `start_time`

  :   The first time to be plotted. Default is the beginning of the time
      series.

  `end_time`

  :   The last time to be plotted. Default is the end of the time
      series.

  `facet`

  :   A boolean value indicating whether to facet the result plot by
      species group. Defaults to TRUE.

  `total`

  :   A boolean value that determines whether the total productivity
      from all species is plotted as well. Default is TRUE.

  `include_repro`

  :   A boolean value that indicates whether to include energy for
      reproduction in productivity estimates. Defaults to using
      [`getEGrowthTime()`](https://cmbeese.github.io/mizerReef/reference/getEGrowthTime.md)
      if FALSE or mizer's
      [`mizer::getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.html)
      if TRUE.

## See also

[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)

Other plotting functions:
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotDegScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegScale.md),
[`plotDegradationScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegradationScale.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md),
[`plotRefugeDensity()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeDensity.md),
[`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md),
[`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md),
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
