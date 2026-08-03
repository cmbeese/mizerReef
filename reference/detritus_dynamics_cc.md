# Detritus dynamics with carrying capacity

Calculates the detritus biomass at the next time step from the current
detritus biomass

## Usage

``` r
detritus_dynamics_cc(params, n, n_other, rates, dt, ...)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
  object

- n:

  A matrix of current species abundances (species x size)

- n_other:

  Other dynamic components.

- rates:

  A list of rates as returned by
  [`getRates()`](https://sizespectrum.org/mizer/reference/getRates.html)

- dt:

  Time step size

- ...:

  Unused

## Value

A single number giving the detritus biomass at next time step

## Details

The time evolution of the detritus biomass \\B\\ is described by

\$\$ \frac{dB_D}{dt} = P_D\left( 1 - \frac{B_D}{K_D} \right) - c_D \\
B_D \$\$

where \\K_D\\ is the system's carrying capacity for detritus in grams/
year, \\c_D\\ is the mass-specific rate of consumption calculated with
[`detritus_consumption()`](https://cmbeese.github.io/mizerReef/reference/detritus_consumption.md)
and \\P_D\\ is the rate at which detritus is produced calculated with
[`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md).

The dynamical equation is solved analytically to

\$\$B_D(t + dt) = B_D(t) \cdot e^{-\frac{dt}{K_D}(P_D+ K_D \\ c_D)} +
\frac{K_D \\ P_D}{P_D+ K_D \\ c_D} \left(1- e^{-\frac{dt}{K_D}(P_D+ K_D
\\ c_D)}\right) \$\$

## See also

[`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md),
[`detritus_consumption()`](https://cmbeese.github.io/mizerReef/reference/detritus_consumption.md),
[`getDetritusConsumption()`](https://cmbeese.github.io/mizerReef/reference/getDetritusConsumption.md),
[`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md)
