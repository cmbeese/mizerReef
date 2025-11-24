# Calculate fisheries productivity for each species group

This function calculates the total potential fisheries productivity by
species group in a given size range using species abundances and the
[`getEGrowthTime()`](https://cmbeese.github.io/mizerReef/reference/getEGrowthTime.md)
function.

## Usage

``` r
getProductivity(
  object,
  time_range = NULL,
  include_repro = FALSE,
  min_fishing_l = NULL,
  max_fishing_l = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `MizerParams` or `MizerSim`. If given a MizerSim
  object, uses the growth rates at the final time of a simulation to
  calculate productivity. If given a MizerParams object, uses the
  initial growth rates.

- time_range:

  The time range (either a vector of values, a vector of min and max
  time, or a single value) to provide productivity for. Default is the
  final time step. Ignored when called with a MizerParams object.

- include_repro:

  A boolean value that indicates whether to include energy for
  reproduction in productivity estimates. Defaults to using
  [`getEGrowthTime()`](https://cmbeese.github.io/mizerReef/reference/getEGrowthTime.md)
  if FALSE or mizer's
  [`mizer::getEReproAndGrowth()`](https://sizespectrum.org/mizer/reference/getEReproAndGrowth.html)
  if TRUE.

- min_fishing_l:

  The minimum length (cm) of fished individuals for productivity
  estimates. Defaults to 7 cm.

- max_fishing_l:

  The maximum length (cm) of fished individuals for productivity
  estimates. Defaults to max length.

- ...:

  Arguments passed on to
  [`mizer::get_size_range_array`](https://sizespectrum.org/mizer/reference/get_size_range_array.html)

  `params`

  :   MizerParams object

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

If called with a MizerParams object, a vector with the productivity in
\\\frac{grams}{m^2}/year\\ for each species group in the model. If
called with a MizerSim object, an array (time x species) containing the
total productivity at each time step for each species.

## Potential fisheries productivity

Productivity refers to the rate at which fish biomass is produced and
available for harvest in a given area over a given period of time.
Productivity cannot be measured in situ.

The productivity \\P_i(w)\\ of species group \\i\\ is given by

\$\$P_i(w) = \int_w^{w+dw} \left( N_i(w) + g_i(w) \right) w \\ dw.\$\$

\\N_i(w)\\ is the abundance density \\(no./m^{2})\\ and \\g_i(w)\\ is
the energy rate available for growth after metabolism and movement have
been accounted for \\(grams/year)\\. The productivity is calculated for
all fish in the size range between `min_fishing_length` and
`max_fishing_length`. These lengths can be the same for all groups or
can be specified as a vector with one value for each species in the
model. The minimum length defaults to \\7 cm\\ regardless of species
group and maximum length defaults to the maximum weight in the model.

## See also

[`getEGrowthTime()`](https://cmbeese.github.io/mizerReef/reference/getEGrowthTime.md),[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)
