# Switch to unstructured resource dynamics with carrying capacities

This function calculates new carrying capacities for algae and detritus
resources and switches the resource dynamics from the standard linear
dynamics to fluxes with a carrying capacity. The carrying capacity is
set at two times the current steady state biomass of the resource.

## Usage

``` r
setURcapacity(params, cap = 1.5, ...)
```

## Arguments

- params:

  a
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
  object

- cap:

  Numeric. Value to scale the steady state biomass by. Defaults to 1.5,
  setting the carrying capacity 50% higher than the current standing
  biomass.

- ...:

  Unused

## Value

A mizer object with updated degradation parameters profiles

## Details

See
[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md)
and
[`detritus_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics_cc.md)
for additional detail.

## See also

[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md),[`detritus_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics_cc.md),
[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
