# Get algae consumption rates

This function returns a named vector with one component for each species
giving the rate in grams/year at which that species consumes algae.

## Usage

``` r
getAlgaeConsumption(params)
```

## Arguments

- params:

  MizerParams

## Value

A named vector with the consumption rates from herbivores

## Details

Unlike
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)
– the mass-specific rate used in
[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md)/[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md)
and in
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md),
which deliberately ignores feeding level, because algal depletion is
modelled as driven by continuous grazing rather than by any individual
consumer's satiation state – this function reports the total,
feeding-level-adjusted rate at which each species actually ingests algae
(e.g. for diagnostic use in
[`plotAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/plotAlgaeConsumption.md)).
For a species with feeding level \\f_i(w)\\, the rate is:

\$\$c\_{i.A} = \rho\_{i.A}\int
w^{m\_{alg}}\\(1-f_i(w))\\N_i(w)\\dw\\B_A\$\$

`f_i(w)` is 0 for species with `satiation = FALSE` (unlimited intake –
by default the carnivores) and follows the usual Holling type II
response otherwise (by default the herbivores/detritivores; see the
`satiation` parameter of
[`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)).
It is a whole-diet quantity computed from a species' total encounter
across all food sources (prey, algae, detritus, plankton combined), not
something computed separately per resource.

## See also

[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md),
[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md),
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)
