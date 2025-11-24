# Detritus dynamics

Calculates the detritus biomass at the next time step based on the
current detritus biomass.

## Usage

``` r
detritus_dynamics(params, n, n_other, rates, dt, ...)
```

## Arguments

- params:

  A MizerParams object

- n:

  A matrix of current species abundances (species x size)

- n_other:

  Other dynamic components.

- rates:

  A list of rates as returned by `getRates()`

- dt:

  Time step size

- ...:

  Unused

## Value

A vector giving the detritus spectrum at the next time step.

## Details

The time evolution of the detritus biomass \\B\\ is described by

\$\$dB_D/dt = P_D - c_D \cdot B_D \$\$

where \\c_D\\ is the mass-specific rate of consumption, calculated with
[`detritus_consumption()`](https://cmbeese.github.io/mizerReef/reference/detritus_consumption.md)and
\\P_D\\ is the rate at which the rest of the system produces detritus
biomass, calculate with
[`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md).

The dynamical equation is solved analytically to

\$\$B_D(t+dt) = B(t)\exp(-c_D \cdot dt) + \frac{P_D}{c_D} (1-\exp(-c_D
\cdot dt)).\$\$

This avoids the stability problems that would arise if we used the Euler
method to solve the equation numerically.

## See also

[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md),
[`detritus_consumption()`](https://cmbeese.github.io/mizerReef/reference/detritus_consumption.md),
[`getDetritusConsumption()`](https://cmbeese.github.io/mizerReef/reference/getDetritusConsumption.md),
[`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md)
