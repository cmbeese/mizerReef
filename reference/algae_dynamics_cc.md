# Algae dynamics with carrying capacity

Calculates the algae biomass at the next time step from the current
algae biomass

## Usage

``` r
algae_dynamics_cc(params, n, n_other, rates, dt, t = 0, ...)
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

- t:

  The current time, used to determine the post-bleaching algae growth
  and capacity boosts, if any (see
  [`getAlgaeBoost()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeBoost.md)).
  Defaults to `0`, which is always before any bleaching event.

- ...:

  Unused

## Value

A single number giving the algae biomass at next time step

## Details

The time evolution of the algae biomass \\B_A(t)\\ is described by

\$\$ \frac{dB_A}{dt} = P_A\left( 1 - \frac{B_A}{K_A} \right) - c_A \\
B_A \$\$

where \\K_A\\ is the system's carrying capacity for algae in grams/
year, \\c_A\\ is the mass-specific rate of consumption calculated with
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)
and \\P_A\\ is the rate at which algae grows, calculated with
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md).

The dynamical equation is solved analytically to

\$\$B_A(t + dt) = B_A(t) \cdot e^{-\frac{dt}{K_A}(P_A+ K_A \\ c_A)} +
\frac{K_A \\ P_A}{P_A+ K_A \\ c_A} \left(1- e^{-\frac{dt}{K_A}(P_A+ K_A
\\ c_A)}\right) \$\$

## See also

[`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md),
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md),
[`getAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeConsumption.md),
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md),
[`getAlgaeBoost()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeBoost.md)
