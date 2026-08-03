# Detritus dynamics

Calculates the detritus biomass at the next time step based on the
current detritus biomass.

## Usage

``` r
detritus_dynamics(params, n, n_other, rates, dt, ...)
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

`B_D(t+dt)` above is a convex combination of `B_D(t)` and `P_D/c_D`, and
is therefore guaranteed non-negative only when `P_D >= 0`. `P_D` (from
[`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md))
is `feces + decomp + external`: `feces` and `decomp` scale with current
fish abundance and mortality, but `external` is a fixed constant, tuned
once at fixed abundances by
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
and never updated again. In a live, non-frozen simulation (e.g. a
multi-year
[`mizer::project()`](https://sizespectrum.org/mizer/reference/project.html)
call with fishing), `feces + decomp` can shrink enough that the fixed
`external` term drives total production negative, breaking the
non-negativity guarantee above and eventually producing negative/`NaN`
detritus biomass. Production is therefore floored at zero
(`production <- max(production, 0)`) before it enters the analytic
update above. This has no effect at any steady state tuned by
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
(there, production equals consumption, which is non-negative by
construction) – it only ever engages once a live simulation has moved
away from that steady state, which is exactly the situation where an
unbounded negative production would otherwise be unphysical: `external`
represents net exchange with the surrounding system (see the "flux of
external detritus is negative" warning in
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md),
documented there as detritus flowing off the reef), and the true floor
on a net production/export rate is zero (nothing can be produced from,
or exported past, an empty pool) rather than an unboundedly negative
constant.

## See also

[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md),
[`detritus_consumption()`](https://cmbeese.github.io/mizerReef/reference/detritus_consumption.md),
[`getDetritusConsumption()`](https://cmbeese.github.io/mizerReef/reference/getDetritusConsumption.md),
[`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md)
