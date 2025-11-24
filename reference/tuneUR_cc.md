# Tune unstructured resources with carrying capacities (algae and detritus) to steady state

For models that use unstructured resources with carrying capacities,
this functions sets the production rates of detritus and algae so that
productions equals consumption at steady state.

## Usage

``` r
tuneUR_cc(params, ...)
```

## Arguments

- params:

  A MizerParams object

- ...:

  unused

## Value

An updated MizerParams object

## Details

With a carrying capacity, the time evolution of the algae biomass
\\B_A(t)\\ is described by

\$\$ \frac{dB_A}{dt} = P_A\left( 1 - \frac{B_A}{K_A} \right) - c_A \\
B_A \$\$

where \\K_A\\ is the system's carrying capacity for algae in grams/
year, \\c_A\\ is the mass-specific rate of consumption calculated with
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)
and \\P_A\\ is the rate at which algae grows, calculated with
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md).

In this tuning function, the growth of rate of algae is set to \\(c_A
\cdot B_A)/(1-\frac{B_A}{K_A})\\ grams per meter squared per year so
that consumption is equal to production for steady state.

Similarly, the time evolution of the detritus biomass \\B_D(t)\\ is
described by

\$\$ \frac{dB_D}{dt} = P_D\left( 1 - \frac{B_D}{K_D} \right) - c_D \\
B_D \$\$

where \\K_D\\ is the system's carrying capacity for detritus in grams/
year, \\c_D\\ is the mass-specific rate of consumption calculated with
[`detritus_consumption()`](https://cmbeese.github.io/mizerReef/reference/detritus_consumption.md)
and \\P_D\\ is the rate at which detritus is produced calculated with
[`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md).
Total detritus production is given with

\$\$p_D = p\_{D.f} + p\_{D.d} + p\_{D.ext}\$\$

In this tuning function, the external production of detritus is set to
\\(c_D \cdot B_D)/(1-\frac{B_D}{K_D}) - P\_{D.f} - P\_{D.d}\\ grams per
meter squared per year so that production equals consumption at steady
state.

## See also

[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md),
[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md),
[`detritus_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics_cc.md)
