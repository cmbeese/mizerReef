# Prepare a steady state model for projections with degradation

This function stores degradation parameters in the `refuge_params` slot
of the params object, and algae adjustment parameters in
`other_params$algae`.

## Usage

``` r
setDegradation(
  params,
  trajectory = NULL,
  deg_scale,
  bleach_time = 2,
  degrade = FALSE,
  algae_boost = FALSE,
  algae_growth_boost = c(1.11, 1.11, 1.11, 1.11),
  algae_capacity_boost = c(2),
  ...
)
```

## Arguments

- params:

  a mizer object

- trajectory:

  Optional character string naming the degradation scenario. Used for
  documentation and identification. No functional effect on the
  simulation.

- deg_scale:

  A matrix (refuge bins x time) giving scaling factors for refuge
  density. Column 1 represents the bleaching year (initial impact), and
  subsequent columns represent years 1, 2, 3... post-bleaching. Values
  are multiplied by the previous timestep's refuge density. Default
  scaling matrices for 15 years with "rubble", "algae", and "recovery"
  trajectories are included as data objects in the package.

- bleach_time:

  The year of the simulation to implement bleaching. Defaults to year 2.

- degrade:

  Logical. Whether to enable habitat degradation during projections,
  i.e. whether
  [`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
  applies the `deg_scale` trajectory once `t` reaches `bleach_time`.
  Default FALSE.

- algae_boost:

  Logical. Should algae growth and/or carrying capacity be adjusted in
  response to bleaching? Default is FALSE.

- algae_growth_boost:

  Numeric vector. Multipliers for algae growth rate at bleaching
  (element 1) and during post-bleaching years (subsequent elements).
  Only used if algae_boost = TRUE. Vector length determines duration:
  length 4 means bleaching year plus 3 post-bleaching years. Default is
  c(1.11, 1.11, 1.11, 1.11) for 4 years of boosting.

- algae_capacity_boost:

  Numeric vector. Multipliers for algae carrying capacity at bleaching
  (element 1) and during post-bleaching years (subsequent elements).
  Only used if algae_boost = TRUE. If shorter than algae_growth_boost,
  will be padded with 1s (no change). Default is c(2.0) (only boosts at
  bleaching year).

- ...:

  Unused

## Value

A mizer object with updated degradation parameters stored in
`refuge_params` and algae adjustment parameters in `other_params$algae`.

## See also

[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md),[`detritus_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics_cc.md),
[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md),
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
