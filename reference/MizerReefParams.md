# Deprecated: Use `newReefParams()` instead

`MizerReefParams()` is a deprecated constructor. Use
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
to create a mizerReef params object.

## Usage

``` r
MizerReefParams(
  ...,
  refuge_params = list(),
  algae_params = list(),
  detritus_params = list()
)
```

## Arguments

- ...:

  Arguments passed to
  [`MizerParams()`](https://sizespectrum.org/mizer/reference/MizerParams.html)

- refuge_params:

  Ignored (now stored in `other_params`)

- algae_params:

  Ignored (now stored in `other_params`)

- detritus_params:

  Ignored (now stored in `other_params`)

## Value

A `mizerReef` object

## See also

[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
