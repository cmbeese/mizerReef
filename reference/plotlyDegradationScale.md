# Interactive Plotly version of plotDegradationScale

Returns an interactive heatmap for degradation scaling parameters.

## Usage

``` r
plotlyDegradationScale(
  object = NULL,
  trajectory = NULL,
  return_data = FALSE,
  ...
)
```

## Arguments

- object:

  An optional object of class
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
  or [MizerSim](https://sizespectrum.org/mizer/reference/MizerSim.html).
  If provided, the function will attempt to extract the degradation
  scaling or trajectory from its `@other_params$refuge_params` slot.

- trajectory:

  Optional. Either a character string (`"rubble"`, `"algae"`,
  `"recovery"`) to use built-in data, or a user-provided numeric
  matrix/data.frame with refuge size bins as rows and bleaching years as
  columns (column 1 = bleaching year, remaining columns = years 1, 2,
  3... post-bleaching).

- return_data:

  Logical. If TRUE, returns the formatted data frame instead of the
  plot. Default FALSE.

- ...:

  Extra parameters passed on to
  [`plotly::ggplotly()`](https://rdrr.io/pkg/plotly/man/ggplotly.html).

## Value

A plotly object
