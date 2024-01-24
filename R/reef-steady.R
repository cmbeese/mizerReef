# Overwrite mizer's steady() function to also set the detritus and algae

#' Project a mizerReef model to steady state
#' 
#' This function tunes the detritus and algae parameters after running
#' mizer's default projectToSteady function.
#' 
#' @param params A \linkS4class{MizerParams} object
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
#' @inheritDotParams tuneUR
#' 
#' @return An object of type \linkS4class{MizerParams}
#' @concept setup
#' @include reef-components.R
#' @export
reefSteady <- function(params, d_func = NULL,
                       t_max = 100, t_per = 1.5, dt = 0.1,
                       tol = 0.1 * dt, return_sim = FALSE,
                       preserve = c("reproduction_level", "erepro", "R_max"),
                       progress_bar = TRUE,...) {

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

    # Force other components to stay at current level
    old_other_dynamics <- params@other_dynamics
    for (res in names(params@other_dynamics)) {
        params@other_dynamics[[res]] <- "constant_dynamics"
    }
    
    if(is.null(d_func)){d_func = distanceSSLogN}

    object <- mizer::projectToSteady(params, distance_func = d_func,
                                     t_per = t_per, t_max = t_max,
                                     dt = dt, tol = tol,
                                     return_sim = return_sim,
                                     progress_bar = progress_bar)

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
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other
    rates <- mizer::getRates(params)
    
    # algae and detritus ----
    params <- tuneUR(params = params, ...)

    if (preserve == "reproduction_level") {
        params <- mizer::setBevertonHolt(params,
                                  reproduction_level = old_reproduction_level)
    } else if (preserve == "R_max") {
        params <- mizer::setBevertonHolt(params,
                                  R_max = old_R_max)
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
