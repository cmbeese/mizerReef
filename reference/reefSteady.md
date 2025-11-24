# Project a mizerReef model to steady state

This function tunes the detritus and algae parameters after running
mizer's default projectToSteady function.

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
  ...
)
```

## Arguments

- params:

  A MizerParams object

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
  reproductive efficiency `erepro`. See `setBevertonHolt()` for an
  explanation of the `reproduction_level`.

- progress_bar:

  A shiny progress object to implement a progress bar in a shiny app.
  Default FALSE.

- ...:

  Arguments passed on to
  [`tuneUR`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)

  :   

## Value

An object of type MizerParams
