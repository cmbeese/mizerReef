# Plot the relative difference in between the total fishable biomasses of each each Species Group at steady state

This functions creates a barplot with the relative change in biomass of
each Species Group within a size range between either (1) two different
mizerParams objects (two models) or (2) two different size ranges.

## Usage

``` r
plotTotalBiomassRelative(
  object1,
  object2,
  diff_method,
  species = NULL,
  min_fishing_l1 = NULL,
  max_fishing_l1 = NULL,
  min_fishing_l2 = NULL,
  max_fishing_l2 = NULL,
  return_data = FALSE,
  ...
)

plotlyTotalBiomassRelative(object1, object2, diff_method, ...)
```

## Arguments

- object1:

  First MizerParams or MizerSim object.

- object2:

  Second MizerParams or MizerSim object.

- diff_method:

  The method to calculate the relative change between models. If
  `percent.change`, the percent change is calculated relative to the
  value from object 1 with formula 100\*(new-old)/old. If `rel.diff` the
  relative difference is returned given by (new - old)/(old + new).

- species:

  The groups to be selected. Optional. By default all target groups are
  selected. A vector of groups names, or a numeric vector with the
  groups indices, or a logical vector indicating for each group whether
  it is to be selected (TRUE) or not.

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

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default value is
  FALSE.

- ...:

  Arguments passed on to
  [`plotTotalBiomass`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)

  `object`

  :   An object of class MizerParams

  `min_fishing_l`

  :   The minimum length (cm) for biomass estimates. Defaults to
      smallest size.

  `max_fishing_l`

  :   The maximum length (cm) of for biomass estimates. Defaults to max
      length.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame with the the total steady state biomass for each functional group
by model is returned as well as another column called `rel_diff`that
gives the relative difference between the two values.

## Details

This function is usually used in conjunction with
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md)
to check for decoupling between biomass and productivity.

The individual productivity rates are calculated by the
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)
function which is passed all additional arguments you supply. See
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)
for more details.

To compare between different size ranges, use the `min_fishing_l1` and
`max_fishing_l1` arguments for the first size range and the
`min_fishing_l2`and `max_fishing_l2` arguments for the second.

## See also

[`plotBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotBiomass.md),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
`plotTotalBiomassRelative()`,
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
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
