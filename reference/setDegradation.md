# Prepare a steady state model for projections with degradation

This function stores degradation parameters in the other_params slot of
the params object.

## Usage

``` r
setDegradation(params, trajectory, deg_scale, bleach_time = 2, ...)
```

## Arguments

- params:

  a mizer object

- trajectory:

  The trajectory for degraded reefs. Options are "rubble", "algae", or
  "recovery"

- deg_scale:

  A 2 x 2 array (refuge size x years post bleaching) that gives the
  values for scaling the refuge density at each size bin. Default
  scaling matrices for 15 years with the "rubble", "algae", and
  "recovery" trajectories are included as data objects in the package.

- bleach_time:

  The year of the simulation to implement bleaching. Defaults to year 2.

- ...:

  Unused

## Value

A mizer object with updated degradation parameters profiles

## See also

[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md),[`detritus_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics_cc.md),
`tune_UR_cc()`,
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
