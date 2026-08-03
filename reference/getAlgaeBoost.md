# Post-bleaching boost multiplier for algae growth or capacity

Computes the multiplier that
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md)
and
[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md)
apply to the baseline algae growth rate or carrying capacity to
represent a temporary post-bleaching boost in algal productivity, as
configured by the `algae_boost`, `algae_growth_boost` and
`algae_capacity_boost` arguments to
[`setDegradation()`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md).

## Usage

``` r
getAlgaeBoost(params, t, boost_vector)
```

## Arguments

- params:

  A MizerParams object

- t:

  The current time

- boost_vector:

  The boost vector to use – `algae_growth_boost` or
  `algae_capacity_boost` from `params@other_params$algae`. Element 1
  gives the multiplier for the bleaching year itself, element 2 for one
  year post-bleaching, and so on. If `t` is later than the last element
  of `boost_vector`, the cumulative product stays fixed at the product
  of the whole vector (i.e. the boost plateaus rather than continuing to
  compound).

## Value

A single number: `1` (no boost) if `algae_boost` is not `TRUE`,
`boost_vector` is empty, or `t` is before the bleaching year; otherwise
the cumulative product of `boost_vector` up to the current post-bleach
year.

## Details

This is a pure function of `t`: it returns the cumulative product of
`boost_vector` up to the element for the current year post-bleaching,
applied fresh to the untouched baseline value every time it is called.
This mirrors how
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
recomputes refuge density fresh from the baseline and `deg_scale` at
every timestep, rather than by mutating `params`: because a single
[`mizer::project()`](https://sizespectrum.org/mizer/reference/project.html)
call treats `params` as read-only across the whole simulation (only
`n`/`n_pp`/`n_other` are threaded forward between timesteps), a rate
function cannot rely on a mutation to `params` persisting into the next
timestep.

## See also

[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md),
[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md),
[`setDegradation()`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md),
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
