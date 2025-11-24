# Rescale algae biomass without changing anything else

This multiplies the algae biomass by a factor and divides the
interaction between all species and algae by the same factor, so as to
keep the total consumption of algae unchanged.

## Usage

``` r
rescale_algae(params, factor)
```

## Arguments

- params:

  A MizerParams object

- factor:

  A number

## Value

An updated MizerParams object
