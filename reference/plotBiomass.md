# Plot the biomass of Species Groups and unstructured components through time

After running a projection, the biomass of each Species Group and each
unstructured component can be plotted against time. The biomass is
calculated within user defined size limits (min_w, max_w, min_l, max_l,
see `get_size_range_array()`).

## Usage

``` r
plotBiomass(
  sim,
  species = NULL,
  start_time,
  end_time,
  y_ticks = 6,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  highlight = NULL,
  return_data = FALSE,
  min_length = NULL,
  ...
)

plotlyBiomass(
  sim,
  species = NULL,
  start_time,
  end_time,
  y_ticks = 6,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  highlight = NULL,
  min_length = NULL,
  ...
)
```

## Arguments

- sim:

  An object of class MizerSim

- species:

  The groups to be selected. Optional. By default all target groups are
  selected. A vector of groups names, or a numeric vector with the
  groups indices, or a logical vector indicating for each group whether
  it is to be selected (TRUE) or not.

- start_time:

  The first time to be plotted. Default is the beginning of the time
  series.

- end_time:

  The last time to be plotted. Default is the end of the time series.

- y_ticks:

  The approximate number of ticks desired on the y axis.

- ylim:

  A numeric vector of length two providing lower and upper limits for
  the y axis. Use NA to refer to the existing minimum or maximum. Any
  values below 1e-20 are always cut off.

- total:

  A boolean value that determines whether the total biomass from species
  is plotted as well. Default is FALSE.

- background:

  A boolean value that determines whether background species are
  included. Ignored if the model does not contain background species.
  Default is TRUE.

- highlight:

  Name or vector of names of the species to be highlighted.

- return_data:

  A boolean value that determines whether the formatted data used for
  the plot is returned instead of the plot itself. Default value is
  FALSE

- min_length:

  minimum length of organism to include in biomass calculation

- ...:

  Arguments passed on to
  [`mizer::get_size_range_array`](https://sizespectrum.org/mizer/reference/get_size_range_array.html)

  `min_w`

  :   Smallest weight in size range. Defaults to smallest weight in the
      model.

  `max_w`

  :   Largest weight in size range. Defaults to largest weight in the
      model.

  `min_l`

  :   Smallest length in size range. If supplied, this takes precedence
      over `min_w`.

  `max_l`

  :   Largest length in size range. If supplied, this takes precedence
      over `max_w`.

## Value

A ggplot2 object, unless `return_data = TRUE`, in which case a data
frame with the four variables 'Year', 'Biomass', 'Species', 'Legend' is
returned.

## See also

`plotBiomass()`,
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md)

Other plotting functions:
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md),
[`plotRefuge()`](https://cmbeese.github.io/mizerReef/reference/plotRefuge.md),
[`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md),
[`plotSpectraRelative()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraRelative.md),
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
