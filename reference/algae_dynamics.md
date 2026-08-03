# Algae dynamics

Calculates the algal biomass at the next time step from the current
algae biomass

## Usage

``` r
algae_dynamics(params, n, n_other, rates, dt, t = 0, ...)
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
  boost, if any (see
  [`getAlgaeBoost()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeBoost.md)).
  Defaults to `0`, which is always before any bleaching event.

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

`B_A(t+dt)` above is a convex combination of `B_A(t)` and `P_A/c_A`, and
is therefore guaranteed non-negative only when `P_A >= 0`. `P_A` (from
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md))
is a fixed, literature-informed production rate representing real algal
primary production, which is not driven by consumer demand: it is set
once (by
[`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)/[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md))
and does not respond to the changing fish abundances of a live,
non-frozen simulation, or get retuned to match consumption – only the
algae *biomass* moves, converging on `P_A/c_A` as `c_A` (and therefore
grazing pressure) changes.
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
tune the algae *biomass* to this same fixed point for a chosen set of
abundances, rather than tuning `P_A`. Production is floored at zero
(`production <- max(production, 0)`) before it enters the analytic
update above, so that this non-negativity guarantee holds even if some
future/degraded parameterisation ever drove `P_A` negative. This has no
effect at any steady state tuned by
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
(there, `P_A` is a literature-informed constant and is non-negative by
construction) – it only ever engages away from that steady state.

## See also

[`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md),
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md),
[`getAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeConsumption.md),
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md),
[`getAlgaeBoost()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeBoost.md)
