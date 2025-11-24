# Get algae consumption rates

This function returns a named vector with one component for each species
giving the rate in grams/year at which that species consumes algae

## Usage

``` r
getAlgaeConsumption(params)
```

## Arguments

- params:

  MizerParams

## Value

A named vector with the consumption rates from herbivores

## Algae consumption

The rate at which herbivorous consumer groups encounter algae biomass
\\E\_{i.A}(w)\\ is controlled by the parameter \\\rho\_{A.i}\\. It
scales with the size of the consumer raised to an allometric exponent
\\m\_{alg}\\ which is taken from empirical data.

\$\$E\_{i.A}(w)=\rho\_{i.A}\\ w^{m\_{alg}}\\B_A\$\$

The mass specific consumption rate then accounts for the preference of
functional group \$i\$ for algae, \\\theta\_{i.A}\\. This gives the
mass-specific algae consumption rate:

\$\$c_A = \sum_i\int\rho\_{i.A}\\ w^{m\_{alg}}
N_i(w)\theta\_{i.A}\\dw\$\$

## See also

[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md),
[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md),
`getAlgaeConsumption()`
