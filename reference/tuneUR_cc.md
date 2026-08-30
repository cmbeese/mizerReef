# Tune unstructured resources with carrying capacities (algae and detritus) to steady state

For models that use unstructured resources with carrying capacities,
this function sets the algae biomass and the production rate of detritus
so that each is at steady state.

## Usage

``` r
tuneUR_cc(params, info_level = mizer::default_info_level(), ...)
```

## Arguments

- params:

  A MizerParams object

- info_level:

  How much mizer should say about the choices it makes here. Level 1
  keeps only the reports that tell you something went differently from
  how you asked; 0 is silence. See
  [`mizer::default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.html).

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

Unlike detritus, algal production on a reef is primary production and is
not driven by consumer demand: it is a property of the algae itself (see
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md)),
left unchanged by this function. Instead, this tuning function solves
the above steady-state condition for the algae biomass, setting it to
\\B_A = (K_A \cdot P_A)/(P_A + K_A \cdot c_A)\\ grams per meter squared,
so that production equals consumption for the given, fixed production
rate. This means that, all else equal, a *decrease* in grazing pressure
(lower \\c_A\\) increases the tuned algae biomass rather than reducing
\\P_A\\ to compensate.

The time evolution of the detritus biomass \\B_D(t)\\ is described by

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
state. Unlike algae, detritus production genuinely is driven by the rest
of the system (fish egestion, senescence) plus an external flux, so
tuning the (external) production rate to match a chosen biomass remains
appropriate.

## See also

[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md),
[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md),
[`detritus_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics_cc.md)
