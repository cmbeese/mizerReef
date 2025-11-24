# Tune unstructured resources (algae and detritus) to steady state

This first sets the rate of degradation of algae so that for the given
abundances, the algae is at steady state. It then sets the rate at which
detritus flows in from external sources (e.g. the pelagic zone) so that
for the given abundances the detritus is at steady state.

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

## See also

[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md),
[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md),
[`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md)
