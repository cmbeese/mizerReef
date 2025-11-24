# Rescale algae and detritus biomass without changing anything else

This multiplies the algae & detritus biomass by a factor and divides the
interaction between all groups and unstructured resource by the same
factor, so as to keep the total consumption of these resources
unchanged.

## Usage

``` r
rescaleComponents(params, algae_factor = 1, detritus_factor = 1)
```

## Arguments

- params:

  A MizerParams object

- algae_factor:

  A number to scale algae by

- detritus_factor:

  A number to scale detritus by

## Value

An updated MizerParams object
