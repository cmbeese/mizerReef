# Mass-specific detritus consumption rate

This mass-specific consumption rate is used in
[`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md)
to calculate the detritus biomass at the next time step. To get the
non-mass-specific consumption rate, use
[`getDetritusConsumption()`](https://cmbeese.github.io/mizerReef/reference/getDetritusConsumption.md).

## Usage

``` r
detritus_consumption(params, n = params@initial_n, rates = getRates(params))
```

## Arguments

- params:

  MizerParams

- n:

  A matrix of current species abundances (species x size)

- rates:

  A list of rates as returned by
  [`getRates()`](https://sizespectrum.org/mizer/reference/getRates.html)

## Value

The mass-specific consumption rate of detritus in grams per year.

## Details

The rho parameter for each functional group is stored in
`other_params(params)$detritus$rho`

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

## Examples

``` r
data(caribbean_3_model)
detritus_consumption(caribbean_3_model)
#> [1] 3.146801e+12
```
