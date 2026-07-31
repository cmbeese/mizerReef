# Overwrite mizer's steady() function to also set the detritus and algae

#' Project a mizerReef model to steady state
#'
#' This function tunes the detritus and algae biomass after running mizer's
#' default `projectToSteady()` function on the fish sub-model.
#'
#' @details
#' Algae and detritus are treated differently while the fish sub-model
#' converges. Detritus is frozen at its current biomass throughout (its
#' `other_dynamics` is temporarily replaced with [constant_dynamics()]),
#' because its production genuinely depends on a not-yet-tuned external
#' flux (see [tuneUR()]/[tuneUR_cc()], [getDetritusProduction()]) that
#' isn't set until the very end. Algae, by contrast, is left live and
#' evolving via its own dynamics function
#' ([algae_dynamics()]/[algae_dynamics_cc()]) throughout the loop, because
#' algae production is a fixed, literature-informed constant (see
#' [setAlgaeParams()]) that never needs retuning -- letting it co-adapt
#' with the fish sub-model as it converges lets the two reach a genuine
#' joint steady state together. (This matters because mizer's own
#' convergence check, `distanceSSLogN()`, only looks at fish abundances,
#' not at the resource pools -- freezing algae and then moving it in one
#' shot at the very end, as used to happen, could leave fish not actually
#' adapted to the final algae biomass.) Algae is frozen instead, like
#' detritus, when `new_refuge == TRUE`, matching the fact that the final
#' tuning step is also skipped in that case (see `new_refuge` in
#' [newReefParams()]).
#'
#' After the fish sub-model converges, [tuneUR()] or [tuneUR_cc()] (chosen
#' based on whether the model uses the carrying-capacity resource
#' formulation) is called once to bring detritus to its own steady state
#' and to make algae's (already live-evolved) biomass exact, unless
#' `new_refuge == TRUE`, in which case algae and detritus are left
#' completely untouched.
#'
#' @param params A [MizerParams] object
#'
#' @param d_func    Optional. A function that will be called after every t_per
#'                  years with both the previous and the new state and that
#'                  should return a number that in some sense measures the
#'                  distance between the states. By default this uses the
#'                  function distanceSSLogN() that you can use as a model
#'                  for your own distance function.
#'
#' @param t_max The maximum number of years to run the simulation. Default is 100.
#'
#' @param t_per The simulation is broken up into shorter runs of `t_per` years,
#'              after each of which we check for convergence. Default value is
#'              1.5. This should be chosen as an odd multiple of the timestep
#'              `dt` in order to be able to detect period 2 cycles.
#'
#' @param dt The time step to use in `project()`.
#'
#' @param tol   The simulation stops when the relative change in the egg
#'              production RDI over `t_per` years is less than `tol` for every
#'              species.
#'
#' @param return_sim    If TRUE, the function returns the MizerSim object
#'                      holding the result of the simulation run. If FALSE
#'                      (default) the function returns a MizerParams object
#'                      with the "initial" slots set to the steady state.
#'
#' @param preserve `r lifecycle::badge("experimental")`
#'   Specifies whether the `reproduction_level` should be preserved (default)
#'   or the maximum reproduction rate `R_max` or the reproductive
#'   efficiency `erepro`. See [setBevertonHolt()] for an explanation
#'   of the `reproduction_level`.
#'
#' @param progress_bar  A shiny progress object to implement a progress bar in a
#'                      shiny app. Default FALSE.
#'
#' @param ... unused
#'
#' @return An object of type [MizerParams]
#' @concept setup
#' @include reef-components.R
#' @export
reefSteady <- function(params, d_func = NULL,
                       t_max = 100, t_per = 1.5, dt = 0.1,
                       tol = 0.1 * dt, return_sim = FALSE,
                       preserve = c("reproduction_level", "erepro", "R_max"),
                       progress_bar = TRUE, ...) {
    # Check if params are valid
    params <- mizer::validParams(params)

    # Choose values to preserve from old models, can be reproduction level,
    # erepro, or R_max
    preserve <- match.arg(preserve)
    old_reproduction_level <- mizer::getReproductionLevel(params)
    old_R_max <- params@species_params$R_max
    old_erepro <- params@species_params$erepro

    # Force the reproduction to stay at the current level
    params@species_params$constant_reproduction <- getRDD(params)
    old_rdd_fun <- params@rates_funcs$RDD
    params@rates_funcs$RDD <- "constantRDD"

    # Force other components to stay at current level, EXCEPT algae: algae
    # production is a fixed, literature-informed constant that is never
    # retuned to match consumption (see tuneUR()/tuneUR_cc()), and its
    # dynamics (algae_dynamics()/algae_dynamics_cc()) converge to their
    # steady state very fast relative to the fish sub-model, so leaving
    # algae live lets it co-adapt with fish throughout convergence instead
    # of being frozen and then jumped in one shot by the final
    # tuneUR()/tuneUR_cc() call below. That one-shot jump matters because
    # mizer's own convergence check (distanceSSLogN()) only looks at fish
    # abundances, not at n_other -- so a frozen-then-jumped algae biomass
    # left fish never having had a chance to adapt to the jump, and the
    # "steady state" was not actually a joint fixed point. Only skip
    # freezing algae when the final tuning step will actually run (i.e.
    # when new_refuge == FALSE); when new_refuge == TRUE, algae is meant to
    # be left completely untouched (see newReefParams()'s new_refuge docs).
    old_other_dynamics <- params@other_dynamics
    for (res in names(params@other_dynamics)) {
        if (res == "algae" && params@other_params$new_refuge == FALSE) {
            next
        }
        params@other_dynamics[[res]] <- "constant_dynamics"
    }

    if (is.null(d_func)) {
        d_func <- distanceSSLogN
    }

    object <- mizer::projectToSteady(params,
        distance_func = d_func,
        t_per = t_per, t_max = t_max,
        dt = dt, tol = tol,
        return_sim = return_sim,
        progress_bar = progress_bar
    )

    if (return_sim) {
        params <- object@params
    } else {
        params <- object
    }

    # Restore original RDD and other dynamics
    params@rates_funcs$RDD <- old_rdd_fun
    params@other_dynamics <- old_other_dynamics
    params@species_params$constant_reproduction <- NULL

    # bring algae and detritus back into steady state
    # n <- params@initial_n
    # n_pp <- params@initial_n_pp
    # n_other <- params@initial_n_other
    # rates <- mizer::getRates(params)

    # algae and detritus ----
    if (params@other_params$new_refuge == FALSE) {
        # use_UR_cc is only set once setAlgaeParams()/setDetritusParams()
        # have been called with use_UR_cc = TRUE; it is unset (NULL) on
        # models that have never opted into the carrying-capacity-scaled
        # formulation, so treat that as FALSE rather than erroring.
        cc <- isTRUE(params@other_params$use_UR_cc)
        if (cc) {
            params <- tuneUR_cc(params = params, ...)
        } else {
            params <- tuneUR(params = params, ...)
        }
    }

    if (preserve == "reproduction_level") {
        params <- mizer::setBevertonHolt(params,
            reproduction_level = old_reproduction_level
        )
    } else if (preserve == "R_max") {
        params <- mizer::setBevertonHolt(params,
            R_max = old_R_max
        )
    } else {
        params <- mizer::setBevertonHolt(params, erepro = old_erepro)
    }

    if (return_sim) {
        object@params <- params
        return(object)
    } else {
        params@time_modified <- lubridate::now()
        return(params)
    }
}

environment(reefSteady) <- asNamespace("mizer")
utils::assignInNamespace("steady", reefSteady, ns = "mizer")
