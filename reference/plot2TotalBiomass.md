# Plot the total biomass of two models or of two different size ranges in the same plot

When called with a
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
object the steady state biomasses are plotted.

## Usage

``` r
plot2TotalBiomass(
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

plotly2TotalBiomass(object1, object2, species = NULL, ...)
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
  [`plotTotalBiomass`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)

  `object`

  :   An object of class
      [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
      or
      [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html).
      If a
      [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html)
      object is provided, the biomass at the last time step is used.

  `min_fishing_l`

  :   The minimum length (cm) for biomass estimates. Defaults to
      smallest size.

  `max_fishing_l`

  :   The maximum length (cm) of for biomass estimates. Defaults to max
      length.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame with the the total steady state biomass for each species by model
is returned as well as another column called `rel_diff`that gives the
relative difference between the two values.

## See also

[`plotBiomass()`](https://sizespectrum.org/mizer/reference/plotBiomass.html),
`plot2TotalBiomass()`,
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md)

Other plotting functions:
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
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
