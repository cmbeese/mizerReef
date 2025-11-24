# Stepped refuge profile for tuning steady states

This is a 2-dimensional array containing start and end lengths for size
bins and prop_protect decreasing from 30% to 10% over the ten bins.

## Usage

``` r
step_tune
```

## Format

dataframe

## Source

PhD Thesis

## Details

This profile provides more protection to smaller size classes than
larger ones, as would be observed on a natural reef.

These refuge parameters are intended for tuning the steady state when
using the density-dependent competitive method. The tuning profile
provides a constant proportion of refuges to all fish up to 50 cm in
length.

When creating a model using the competitive method, you should run
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
with the "binned" method and thing tuning profile.

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
