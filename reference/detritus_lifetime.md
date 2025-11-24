# Expected detritus lifetime

The expected detritus lifetime is defined as the inverse of the
mass-specific detritus consumption rate.

## Usage

``` r
detritus_lifetime(params)

detritus_lifetime(params) <- value
```

## Arguments

- params:

  A MizerParams object

- value:

  A number with the new value for the expected lifetime in years

  Assigning a new value to the detritus lifetime rescales the detritus
  abundance while keeping the total consumption of detritus the same (by
  adjusting the interaction strength of species with detritus).

## Value

The number giving the expected lifetime in years.
