# Constant refuge profile for tuning steady states

This is a 2-dimensional array containing start and end lengths for size
bins and `prop_protect` equal to 20% for all size bins up to 50 cm in
length.

## Usage

``` r
tuning_profile
```

## Format

dataframe

## Source

PhD Thesis

## Details

These refuge parameters are intended for tuning the steady state when
using the density-dependent competitive method. The tuning profile
provides a constant proportion of refuges to all fish up to 50 cm in
length.

When creating a model using the competitive method, you should run
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
with the "binned" method and a proportional tuning profile.

After species biomasses and growth rates have been calibrated to match
empirical observations, use the
[`newRefuge()`](https://cmbeese.github.io/mizerReef/reference/newRefuge.md)
function to implement your competitive refuge parameters. After using
[`newRefuge()`](https://cmbeese.github.io/mizerReef/reference/newRefuge.md),
make sure to iterate through
[`mizer:: matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.html),
[`matchReefGrowth()`](https://cmbeese.github.io/mizerReef/reference/matchReefGrowth.md),
and
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
again to regain the steady state.
