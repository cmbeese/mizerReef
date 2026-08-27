# Reef-aware steady state: mizer's steady-state machinery on the fish
# sub-model, plus the algae and detritus tuning that mizer knows nothing
# about.

#' Project a mizerReef model to steady state
#'
#' This function tunes the detritus and algae biomass after running mizer's
#' [mizer::findSteadyState()] on the fish sub-model.
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
#' `reefSteady()` is also registered as the [mizer::steady()] and
#' [mizer::tuneSteadyState()] method for `mizerReef` objects, so calling
#' either of those on a reef model does the reef-aware thing. Earlier
#' versions instead replaced `mizer::steady()` in mizer's namespace, which
#' broke `steady()` for every non-reef model in the session.
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
#' @param ... Passed on to [mizer::findSteadyState()] (or, when
#'   `return_sim = TRUE`, to [mizer::projectUntilSettled()]), so that
#'   arguments such as `effort`, `method` or `info_level` can be given.
#'
#' @param info_level How much mizer should say about the choices it makes
#'   here. Level 1 keeps only the reports that tell you something went
#'   differently from how you asked; 0 is silence. See
#'   [mizer::default_info_level()].
#' @return An object of type [MizerParams]
#' @concept setup
#' @examples
#' data(caribbean_3_model)
#' params <- reefSteady(caribbean_3_model)
#' @include reef-components.R
#' @export
reefSteady <- function(params, d_func = NULL,
                       t_max = 100, t_per = 1.5, dt = 0.1,
                       tol = 0.1 * dt, return_sim = FALSE,
                       preserve = c("reproduction_level", "erepro", "R_max"),
                       progress_bar = TRUE,
                       info_level = mizer::default_info_level(), ...) {
    # One handler for the whole call, so that the reports raised by
    # tuneUR()/tuneUR_cc() below arrive together with anything mizer's
    # steady-state machinery has to say.
    mizer::with_info_level(info_level = info_level, {
    # Check if params are valid
    params <- mizer::validParams(params)

    # Choose values to preserve from old models, can be reproduction level,
    # erepro, or R_max
    preserve <- match.arg(preserve)
    old_reproduction_level <- mizer::getReproductionLevel(params)
    old_R_max <- params@species_params$R_max
    old_erepro <- params@species_params$erepro

    # Force the reproduction to stay at the current level.
    # `constant_reproduction` is written into the slot on purpose, and taken
    # out again below: it is a temporary flag read by `constantRDD()` for the
    # duration of this call only, exactly as mizer's own
    # `tune_steady_project()` sets it. Routing it through
    # `species_params<-()` would record it as a given species parameter and
    # then withdraw the column again on every call.
    params@species_params$constant_reproduction <- mizer::getRDD(params)
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
        d_func <- mizer::distanceSSLogN
    }

    # mizer 3.3 renamed projectToSteady() to findSteadyState(), with
    # projectUntilSettled() for the trajectory, and renamed t_per/tol to
    # t_check/distance_tol. The new finders stop only once the biomass
    # drift has settled as well (require_steady = TRUE); reefSteady()'s
    # documented stopping rule is the distance function alone, so we ask
    # for the old rule explicitly, exactly as mizer's own steady() and
    # projectToSteady() wrappers do. t_save = t_per keeps the trajectory
    # spacing these functions have always returned.
    #
    # These defaults are merged with `...` via modifyList() (user values
    # win) rather than passed alongside `...`, so a caller can override any
    # of them -- e.g. reefSteady(params, distance_tol = 0.05) -- instead of
    # hitting "formal argument matched by multiple actual arguments" from
    # naming the same argument twice.
    steady_args <- utils::modifyList(
        list(
            distance_func = d_func,
            t_check = t_per, t_max = t_max,
            dt = dt, t_save = t_per, distance_tol = tol,
            require_steady = FALSE,
            progress_bar = progress_bar, info_level = info_level
        ),
        list(...)
    )
    if (return_sim) {
        object <- do.call(mizer::projectUntilSettled,
                           c(list(params), steady_args))
        params <- object@params
    } else {
        object <- do.call(
            mizer::findSteadyState,
            utils::modifyList(list(params = params, solver = "project"),
                               steady_args)
        )
        params <- object
    }

    # Restore original RDD and other dynamics
    params@rates_funcs$RDD <- old_rdd_fun
    params@other_dynamics <- old_other_dynamics
    params@species_params$constant_reproduction <- NULL

    # algae and detritus ----
    if (params@other_params$new_refuge == FALSE) {
        # use_UR_cc is only set once setAlgaeParams()/setDetritusParams()
        # have been called with use_UR_cc = TRUE; it is unset (NULL) on
        # models that have never opted into the carrying-capacity-scaled
        # formulation, so treat that as FALSE rather than erroring.
        cc <- isTRUE(params@other_params$use_UR_cc)
        if (cc) {
            params <- tuneUR_cc(params = params)
        } else {
            params <- tuneUR(params = params)
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
        object
    } else {
        params@time_modified <- lubridate::now()
        params
    }
    })
}

#' Steady state methods for mizerReef models
#'
#' `mizerReef` methods for mizer's steady-state generics, both of which run
#' [reefSteady()] so that the algae and detritus pools are tuned along with
#' the fish sub-model. `mizer::steady()` is the superseded name;
#' `mizer::tuneSteadyState()` is the current one.
#'
#' `tuneSteadyState()` supports only `solver = "project"` here, because
#' `reefSteady()` needs the projection in order to let algae co-adapt with
#' the fish spectra (see [reefSteady()]).
#'
#' @param params A `mizerReef` object
#' @param t_max,t_per,dt,t_save,tol Control the projection. `t_save` is
#'   accepted for compatibility with the generic and, as in mizer's own
#'   method, does not affect the result.
#' @param amplitude_tol,amp_rel_tol,extinction_threshold,method,info_level
#'   Passed on to mizer's steady-state machinery.
#' @param return_sim If TRUE, return the `MizerSim` of the run rather than
#'   the tuned `MizerParams`.
#' @param preserve See [reefSteady()].
#' @param progress_bar A shiny progress object, or FALSE.
#' @param solver Only `"project"` is supported for reef models.
#' @param effort The fishing effort to hold during the run.
#' @param ... Passed on to [reefSteady()].
#'
#' @return A [MizerParams] object, or a `MizerSim` when `return_sim = TRUE`.
#' @seealso [reefSteady()]
#' @concept setup
#' @name reef-steady-methods
NULL

#' @rdname reef-steady-methods
#' @method steady mizerReef
#' @export
steady.mizerReef <- function(params, t_max = 100, t_per = 1.5, dt = 0.1,
                             t_save = dt, tol = 0.1 * dt,
                             amplitude_tol = 0.01, amp_rel_tol = 0.01,
                             extinction_threshold = 1e-6,
                             return_sim = FALSE,
                             preserve = c("reproduction_level", "erepro",
                                          "R_max"),
                             progress_bar = TRUE,
                             info_level = mizer::default_info_level(),
                             method = c("euler", "predictor_corrector",
                                        "tr_bdf2")) {
    # Forward the cycle-detection and integration arguments only when the
    # caller actually gave them. Their defaults on mizer's `steady()` generic
    # are not all the same as the ones `findSteadyState()` uses -- `steady()`
    # documents `amp_rel_tol = 0.01`, the finder uses 0.1 -- and passing the
    # generic's default through would quietly change what `steady()` does to
    # a reef model relative to calling `reefSteady()` directly.
    extra <- list()
    if (!missing(amplitude_tol)) extra$amplitude_tol <- amplitude_tol
    if (!missing(amp_rel_tol)) extra$amp_rel_tol <- amp_rel_tol
    if (!missing(extinction_threshold)) {
        extra$extinction_threshold <- extinction_threshold
    }
    if (!missing(method)) extra$method <- method
    if (!missing(info_level)) extra$info_level <- info_level

    do.call(reefSteady, c(
        list(params,
             t_max = t_max, t_per = t_per, dt = dt, tol = tol,
             return_sim = return_sim, preserve = preserve,
             progress_bar = progress_bar),
        extra
    ))
}

#' @rdname reef-steady-methods
#' @method tuneSteadyState mizerReef
#' @export
tuneSteadyState.mizerReef <- function(params,
                                      solver = c("project", "newton"),
                                      effort = params@initial_effort,
                                      preserve = c("reproduction_level",
                                                   "erepro", "R_max"),
                                      info_level = mizer::default_info_level(),
                                      ...) {
    solver <- match.arg(solver)
    if (solver == "newton") {
        stop("`tuneSteadyState()` on a mizerReef model supports only ",
             "`solver = \"project\"`: reefSteady() lets algae co-adapt with ",
             "the fish spectra as they settle, which the Newton solver ",
             "cannot do. See `?reefSteady`.")
    }
    reefSteady(params,
        preserve = preserve, effort = effort, info_level = info_level, ...
    )
}
