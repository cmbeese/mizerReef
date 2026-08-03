# Calibrate the model scale to match total observed number

Replaces mizer's
[`mizer::calibrateNumber()`](https://sizespectrum.org/mizer/reference/calibrateNumber.html)
function. Given a MizerParams object `params` for which number
observations are available for at least some species via the
`number_observed` column in the species_params data frame, this function
returns an updated MizerParams object which is rescaled with
[`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md)
so that the total number in the model agrees with the total observed
number.

## Usage

``` r
calibrateReefNumber(params)
```

## Arguments

- params:

  A MizerParams object

## Value

A MizerParams object

## Details

Number observations usually only include individuals above a certain
size. This size should be specified in a number_cutoff column of the
species parameter data frame. If this is missing, it is assumed that all
sizes are included in the observed number, i.e., it includes larval
number.

After using this function the total number in the model will match the
total number, summed over all species. However the numbers of the
individual species will not match observations yet, with some species
having numbers that are too high and others too low. So after this
function you may want to use
[`matchNumbers()`](https://sizespectrum.org/mizer/reference/matchNumbers.html).
This is described in the blog post at https://bit.ly/2YqXESV.

If you have observations of the yearly yield instead of numbers, you can
use
[`calibrateYield()`](https://sizespectrum.org/mizer/reference/calibrateYield.html)
instead of this function.
