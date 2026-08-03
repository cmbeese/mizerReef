# Rescale detritus biomass without changing anything else

This multiplies the detritus biomass by a factor and divides the
interaction between all species and the detritus by the same factor, so
as to keep the total consumption of detritus unchanged.

## Usage

``` r
rescale_detritus(params, factor)
```

## Arguments

- params:

  A MizerParams object

- factor:

  A number

## Value

An updated MizerParams object
