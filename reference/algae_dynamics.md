# Algae dynamics

Calculates the algal biomass at the next time step from the current
algae biomass

## Usage

``` r
algae_dynamics(params, n, n_other, rates, dt, ...)
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

A single number giving the algae biomass at next time step

## Details

The time evolution of the algal biomass \\B\\ is described by

\$\$dB_A/dt = P_A - c_A \cdot B_A\$\$

where \\c_A\\ is the mass-specific rate of consumption calculated with
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)
and \\P_A\\ is the rate at which algae grows, calculated with
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md).

The dynamical equation is solved analytically to

\$\$B_A(t+dt) = B_A(t) e^{(- c_A \cdot dt)} +\frac{P_A}{c_A} (1-
e^{(-c_A \cdot dt)}).\$\$

This avoids the stability problems that would arise if we used the Euler
method to solve the equation numerically.

## See also

[`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md),
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md),
[`getAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeConsumption.md),
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md)
