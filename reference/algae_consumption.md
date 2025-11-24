# Mass-specific algae consumption rate

This mass-specific consumption rate is used in
[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md)
to calculate the algae biomass at the next time step. To get the
non-mass-specific consumption rate, use
[`getAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeConsumption.md).

## Usage

``` r
algae_consumption(params, n = params@initial_n, rates = getRates(params))
```

## Arguments

- params:

  MizerParams

- n:

  A matrix of current species abundances (species x size)

- rates:

  A list of rates as returned by `getRates()`

## Value

The mass-specific consumption rate of algae in grams per year.

## Details

The rho parameter for herbivorous fish groups is stored in
`other_params(params)$algae$rho`

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
