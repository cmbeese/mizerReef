# Steady state methods for mizerReef models

`mizerReef` methods for mizer's steady-state generics, both of which run
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
so that the algae and detritus pools are tuned along with the fish
sub-model.
[`mizer::steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
is the superseded name;
[`mizer::tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.html)
is the current one.

## Usage

``` r
# S3 method for class 'mizerReef'
steady(
  params,
  t_max = 100,
  t_per = 1.5,
  dt = 0.1,
  t_save = dt,
  tol = 0.1 * dt,
  amplitude_tol = 0.01,
  amp_rel_tol = 0.01,
  extinction_threshold = 1e-06,
  return_sim = FALSE,
  preserve = c("reproduction_level", "erepro", "R_max"),
  progress_bar = TRUE,
  info_level = mizer::default_info_level(),
  method = c("euler", "predictor_corrector", "tr_bdf2")
)

# S3 method for class 'mizerReef'
tuneSteadyState(
  params,
  solver = c("project", "newton"),
  effort = params@initial_effort,
  preserve = c("reproduction_level", "erepro", "R_max"),
  info_level = mizer::default_info_level(),
  ...
)
```

## Arguments

- params:

  A `mizerReef` object

- t_max, t_per, dt, t_save, tol:

  Control the projection. `t_save` is accepted for compatibility with
  the generic and, as in mizer's own method, does not affect the result.

- amplitude_tol, amp_rel_tol, extinction_threshold, method, info_level:

  Passed on to mizer's steady-state machinery.

- return_sim:

  If TRUE, return the `MizerSim` of the run rather than the tuned
  `MizerParams`.

- preserve:

  See
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md).

- progress_bar:

  A shiny progress object, or FALSE.

- solver:

  Only `"project"` is supported for reef models.

- effort:

  The fishing effort to hold during the run.

- ...:

  Passed on to
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md).

## Value

A
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
object, or a `MizerSim` when `return_sim = TRUE`.

## Details

`tuneSteadyState()` supports only `solver = "project"` here, because
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
needs the projection in order to let algae co-adapt with the fish
spectra (see
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)).

## See also

[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
