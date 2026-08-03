# Tune unstructured resources (algae and detritus) to steady state

This first sets the algae biomass so that, for the given abundances and
the algae's fixed production rate (see
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md)),
algae is at steady state. It then sets the rate at which detritus flows
in from external sources (e.g. the pelagic zone) so that for the given
abundances the detritus is at steady state.

## Usage

``` r
tuneUR(params, ...)
```

## Arguments

- params:

  A MizerParams object

- ...:

  unused

## Value

An updated MizerParams object

## Details

The time evolution of the algae biomass \\B_A(t)\\ is described by
\\dB_A/dt = P_A - c_A \cdot B_A\\, where \\c_A\\ is the mass-specific
consumption rate from
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)
and \\P_A\\ is the algae production rate from
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md).
Unlike detritus, algal production on a reef is primary production, not
driven by consumer demand: \\P_A\\ is a fixed property of the algae,
left unchanged by this function. Instead, this tuning function solves
\\dB_A/dt = 0\\ for the algae biomass, setting it to \\B_A = P_A /
c_A\\. This means that, all else equal, a *decrease* in grazing pressure
(lower \\c_A\\) increases the tuned algae biomass rather than reducing
\\P_A\\ to compensate. If \\c_A = 0\\ (no consumers), no finite
steady-state biomass exists and the algae biomass is left unchanged,
with a warning.

## See also

[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md),
[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md),
[`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md)

## Examples

``` r
data(caribbean_3_model)
params <- tuneUR(caribbean_3_model)
#> Warning: The flux of external detritus is negative.
```
