# Project a mizerReef model to steady state

This function tunes the detritus and algae biomass after running mizer's
[`mizer::findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.html)
on the fish sub-model.

## Usage

``` r
reefSteady(
  params,
  d_func = NULL,
  t_max = 100,
  t_per = 1.5,
  dt = 0.1,
  tol = 0.1 * dt,
  return_sim = FALSE,
  preserve = c("reproduction_level", "erepro", "R_max"),
  progress_bar = TRUE,
  info_level = mizer::default_info_level(),
  ...
)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
  object

- d_func:

  Optional. A function that will be called after every t_per years with
  both the previous and the new state and that should return a number
  that in some sense measures the distance between the states. By
  default this uses the function distanceSSLogN() that you can use as a
  model for your own distance function.

- t_max:

  The maximum number of years to run the simulation. Default is 100.

- t_per:

  The simulation is broken up into shorter runs of `t_per` years, after
  each of which we check for convergence. Default value is 1.5. This
  should be chosen as an odd multiple of the timestep `dt` in order to
  be able to detect period 2 cycles.

- dt:

  The time step to use in `project()`.

- tol:

  The simulation stops when the relative change in the egg production
  RDI over `t_per` years is less than `tol` for every species.

- return_sim:

  If TRUE, the function returns the MizerSim object holding the result
  of the simulation run. If FALSE (default) the function returns a
  MizerParams object with the "initial" slots set to the steady state.

- preserve:

  **\[experimental\]** Specifies whether the `reproduction_level` should
  be preserved (default) or the maximum reproduction rate `R_max` or the
  reproductive efficiency `erepro`. See
  [`setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.html)
  for an explanation of the `reproduction_level`.

- progress_bar:

  A shiny progress object to implement a progress bar in a shiny app.
  Default FALSE.

- info_level:

  How much mizer should say about the choices it makes here. Level 1
  keeps only the reports that tell you something went differently from
  how you asked; 0 is silence. See
  [`mizer::default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.html).

- ...:

  Passed on to
  [`mizer::findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.html)
  (or, when `return_sim = TRUE`, to
  [`mizer::projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.html)),
  so that arguments such as `effort`, `method` or `info_level` can be
  given.

## Value

An object of type
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)

## Details

Algae and detritus are treated differently while the fish sub-model
converges. Detritus is frozen at its current biomass throughout (its
`other_dynamics` is temporarily replaced with
[`constant_dynamics()`](https://cmbeese.github.io/mizerReef/reference/constant_dynamics.md)),
because its production genuinely depends on a not-yet-tuned external
flux (see
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md),
[`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md))
that isn't set until the very end. Algae, by contrast, is left live and
evolving via its own dynamics function
([`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md)/[`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md))
throughout the loop, because algae production is a fixed,
literature-informed constant (see
[`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md))
that never needs retuning – letting it co-adapt with the fish sub-model
as it converges lets the two reach a genuine joint steady state
together. (This matters because mizer's own convergence check,
`distanceSSLogN()`, only looks at fish abundances, not at the resource
pools – freezing algae and then moving it in one shot at the very end,
as used to happen, could leave fish not actually adapted to the final
algae biomass.) Algae is frozen instead, like detritus, when
`new_refuge == TRUE`, matching the fact that the final tuning step is
also skipped in that case (see `new_refuge` in
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)).

After the fish sub-model converges,
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md) or
[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
(chosen based on whether the model uses the carrying-capacity resource
formulation) is called once to bring detritus to its own steady state
and to make algae's (already live-evolved) biomass exact, unless
`new_refuge == TRUE`, in which case algae and detritus are left
completely untouched.

`reefSteady()` is also registered as the
[`mizer::steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
and
[`mizer::tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.html)
method for `mizerReef` objects, so calling either of those on a reef
model does the reef-aware thing. Earlier versions instead replaced
[`mizer::steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
in mizer's namespace, which broke `steady()` for every non-reef model in
the session.

## Examples

``` r
data(caribbean_3_model)
params <- reefSteady(caribbean_3_model)
#> Reached the convergence tolerance after 1.5 years. The biomasses change at up to 1.5e-10 per year.
#> Warning: The flux of external detritus is negative.
```
