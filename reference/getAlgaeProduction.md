# Algae production rate

This is the rate in grams/year/m^-2 at which the system produces algae
biomass. Unlike detritus production, this is real primary production and
is not driven by consumer demand: the baseline rate is a fixed,
literature-informed constant set by
[`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)/[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
(see the `initial_algae_growth` argument there for the default value and
its literature basis) and is left unchanged by
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md),
which instead tune the algae *biomass* to the steady state implied by
this fixed production rate and the current consumption. If
post-bleaching algae boosting has been configured with
[`setDegradation()`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md),
the baseline rate is scaled by the boost factor appropriate for time `t`
(see
[`getAlgaeBoost()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeBoost.md)).

## Usage

``` r
getAlgaeProduction(params, t = 0)
```

## Arguments

- params:

  MizerParams

- t:

  The current time, used to determine the post-bleaching algae growth
  boost, if any (see
  [`getAlgaeBoost()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeBoost.md)).
  Defaults to `0`, which is always before any bleaching event, so no
  boost applies.

## Value

The annual growth rate of algae per square meter

## See also

[`getAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeConsumption.md),
[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md),
[`getAlgaeBoost()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeBoost.md)

## Examples

``` r
data(caribbean_3_model)
getAlgaeProduction(caribbean_3_model)
#> [1] 2000
```
