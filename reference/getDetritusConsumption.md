# Get detritus consumption rates

This function returns a named vector with one component for each species
giving the rate in grams/year at which that species consumes detritus

## Usage

``` r
getDetritusConsumption(params)
```

## Arguments

- params:

  MizerParams

## Value

A named vector with the consumption rates from herbivores

## Detritus consumption

The rate at which detritivorous consumer groups encounter detrital
biomass \\E\_{i.D}(w)\\ is controlled by the parameter \\\rho\_{D.i}\\.
It scales with the size of the consumer raised to an allometric exponent
\\m\_{det}\\ which is taken to be the same as the scaling exponent of
the maximum intake rate for fish consumers.

\$\$E\_{i.D}(w)=\rho\_{i.D}\\ w^{m\_{det}}\\B_D \$\$

The mass specific consumption rate then accounts for the preference of
functional group \$i\$ for detritus, \\\theta\_{i.D}\\ and the feeding
level \\f_i(w)\\. This gives the mass-specific detritus consumption
rate:

\$\$c_D = \sum_i\int\rho\_{i.D}\\ w^{m\_{det}} N_i(w)
\left(1-f_i(w)\right) \theta\_{i.D}\\dw\$\$

## See also

[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md),
[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md),
`getDetritusConsumption()`
