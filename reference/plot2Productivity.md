# Plot the fisheries productivity of two models or two different size ranges in the same plot

When called with a
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
object the steady state productivities are plotted.

## Usage

``` r
plot2Productivity(
  object1,
  object2,
  species = NULL,
  name1 = "First",
  name2 = "Second",
  min_fishing_l1 = NULL,
  max_fishing_l1 = NULL,
  min_fishing_l2 = NULL,
  max_fishing_l2 = NULL,
  stack = FALSE,
  return_data = FALSE,
  ...
)

plotly2Productivity(object1, object2, ...)
```

## Arguments

- object1:

  First MizerParams or MizerSim object.

- object2:

  Second MizerParams or MizerSim object.

- species:

  The species to be selected. Optional. By default all target species
  are selected. A vector of species names, or a numeric vector with the
  species indices, or a logical vector indicating for each group whether
  it is to be selected (TRUE) or not.

- name1:

  An optional string with the name for the first model, to be used in
  the legend. Set to "First" by default.

- name2:

  An optional string with the name for the second model, to be used in
  the legend. Set to "Second" by default.

- min_fishing_l1:

  Optional. The minimum length (cm) of fished individuals for model 2.
  Defaults to 7cm. A parameter passed to
  [`getProductivity()`](https://cmbeese.github.io/mizerReef/reference/getProductivity.md).

- max_fishing_l1:

  Optional. The maximum length (cm) of fished individuals for model 1.
  Defaults to max length. A parameter passed to
  [`getProductivity()`](https://cmbeese.github.io/mizerReef/reference/getProductivity.md).

- min_fishing_l2:

  Optional. The minimum length (cm) of fished individuals for model 2.
  Defaults to 7cm. A parameter passed to
  [`getProductivity()`](https://cmbeese.github.io/mizerReef/reference/getProductivity.md).

- max_fishing_l2:

  Optional. The maximum length (cm) of fished individuals for model 1.
  Defaults to max length. A parameter passed to
  [`getProductivity()`](https://cmbeese.github.io/mizerReef/reference/getProductivity.md).

- stack:

  A boolean value that determines whether bars are separated by species.
  Defaults to FALSE. If true, returns a stacked barplot with the total
  biomass for each group instead of individual bars for each group.
  Useful for comparison between steady states.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default value is
  FALSE.

- ...:

  Arguments passed on to
  [`plotProductivity`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)

  `object`

  :   An object of class
      [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
      or
      [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html)

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

  `include_inverts`

  :   A boolean value that determines whether the "inverts" species
      group is included. Default is FALSE, since invertebrate
      productivity is typically not relevant to fishing yield. Only
      takes effect when `species` is not explicitly provided.

  `include_repro`

  :   A boolean value that indicates whether to include energy for
      reproduction in productivity estimates. Defaults to using
      [`getEGrowthTime()`](https://cmbeese.github.io/mizerReef/reference/getEGrowthTime.md)
      if FALSE or mizer's
      [`mizer::getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.html)
      if TRUE.

  `min_fishing_l`

  :   The minimum length (cm) of fished individuals for productivity
      estimates. Defaults to 7 cm.

  `max_fishing_l`

  :   The maximum length (cm) of fished individuals for productivity
      estimates. Defaults to max length.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame with the the productivity for each Species Group by model is
returned.

## See also

[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.html),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
`plot2Productivity()`,
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md)

Other plotting functions:
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
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
