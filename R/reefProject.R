#' Project mizerReef model forward in time and store refuge profiles
#'
#' Simulates size spectrum dynamics for reef models, storing the refuge profile ("vulnerable") at each time step in a `MizerReefSim` object.
#' This function is identical to [mizer::project()] except it also stores the time-varying refuge profile in the `vulnerable` slot.
#'
#' @param object Either a `MizerParams` or `MizerSim` object.
#' @param effort Fishing effort for each gear through time. See [mizer::project()] for details.
#' @param t_max Number of years to project (default = 100). Ignored if effort is an array.
#' @param dt Time step of the solver (default = 0.1).
#' @param t_save Output storage frequency (default = 1). Ignored if effort is an array.
#' @param t_start Start year of the simulation (default = 0). Ignored if effort is an array or if `object` is a `MizerSim`.
#' @param append If `TRUE`, append new results to previous ones (only relevant if `object` is a `MizerReefSim`). Default = TRUE.
#' @param progress_bar Show a progress bar (logical or shiny Progress object).
#' @param ... Other arguments passed to rate functions.
#'
#' @return A `MizerReefSim` object with abundances and refuge profiles through time.
#' @seealso [mizer::project], [MizerReefSim]
#' @examples
#' # See [mizer::project()] for usage examples.
#' @export
#' @concept simulation
#' @family simulation
#' @importFrom methods new
reefProject <- function(object, effort,
                       t_max = 100, dt = 0.1, t_save = 1, t_start = 0,
                       append = TRUE,
                       progress_bar = TRUE,
                       ...) {
    # Set up initial values (similar to mizer's project)
    if (is(object, "MizerSim")) {
        params <- setInitialValues(object@params, object)
        t_start <- getTimes(object)[idxFinalT(object)]
    } else if (is(object, "MizerParams")) {
        params <- validParams(object)
    } else {
        stop("The `object` argument must be either a MizerParams or a MizerSim object.")
    }
    initial_n <- params@initial_n
    initial_n_pp <- params@initial_n_pp
    initial_n_other <- params@initial_n_other

    times <- seq(t_start, t_start + t_max, by = t_save)
    sim <- MizerReefSim(params, t_dimnames = times)
    sim@n[1, , ] <- initial_n
    sim@n_pp[1, ] <- initial_n_pp
    sim@n_other[1, ] <- unserialize(serialize(initial_n_other, NULL))
    sim@effort <- effort

    # Set up progress bar
    if (inherits(progress_bar, "Progress")) {
        progress_bar$set(message = "Running simulation", value = 0)
        proginc <- 1 / length(times)
    } else if (isTRUE(progress_bar)) {
        pb <- progress::progress_bar$new(
            format = "[:bar] :percent ETA: :eta",
            total = length(times), width = 60)
        pb$tick(0)
    }

    # calculate vulnerablity through time to fin refuge profiles
    sim@vulnerable <- vector("list", length(times))

    # Get functions
    resource_dynamics_fn <- get(params@resource_dynamics)
    other_dynamics_fns <- lapply(params@other_dynamics, get)
    rates_fns <- lapply(params@rates_funcs, get)

    n_list <- list(n = initial_n, n_pp = initial_n_pp, n_other = initial_n_other)
    t <- t_start

    # Loop over time
    for (i in seq(2, length(times))) {
        steps <- round((times[[i]] - t) / dt)
        n_list <- project_simple(
            params, n = n_list$n, n_pp = n_list$n_pp, n_other = n_list$n_other,
            t = t, dt = dt, steps = steps,
            effort = effort[i - 1, ],
            resource_dynamics_fn = resource_dynamics_fn,
            other_dynamics_fns = other_dynamics_fns,
            rates_fns = rates_fns, ...)
        t <- t + steps * dt
        sim@n[i, , ] <- n_list$n
        sim@n_pp[i, ] <- n_list$n_pp
        sim@n_other[i, ] <- unserialize(serialize(n_list$n_other, NULL))
        # Store vulnerable profile for this time step
        sim@vulnerable[[i]] <- reefVulnerable(params, n_list$n, n_list$n_pp, n_list$n_other, t)
        # Advance progress bar
        if (inherits(progress_bar, "Progress")) {
            progress_bar$inc(amount = proginc)
        } else if (isTRUE(progress_bar)) {
            pb$tick()
        }
    }

    # append to previous simulation ----
    if (is(object, "MizerReefSim") && append) {
        no_t_old <- dim(object@n)[1]
        no_t <- length(times)
        new_t_dimnames <- c(as.numeric(dimnames(object@n)[[1]]),
                            times[2:no_t])
        new_sim <- MizerReefSim(params, t_dimnames = new_t_dimnames)
        old_indices <- 1:no_t_old
        new_indices <- seq(from = no_t_old + 1, length.out = no_t - 1)
        new_sim@n[old_indices, , ]  <- object@n
        new_sim@n[new_indices, , ]  <- sim@n[2:no_t, , ]
        new_sim@n_pp[old_indices, ] <- object@n_pp
        new_sim@n_pp[new_indices, ] <- sim@n_pp[2:no_t, ]
        new_sim@n_other[old_indices, ]  <- object@n_other
        new_sim@n_other[new_indices, ]  <- sim@n_other[2:no_t, ]
        new_sim@effort[old_indices, ] <- object@effort
        new_sim@effort[new_indices, ] <- sim@effort[2:no_t, ]
        # Combine vulnerable profiles
        new_sim@vulnerable <- c(object@vulnerable, sim@vulnerable[2:no_t])
        new_sim@vulnerable_initial <- object@vulnerable_initial
        return(new_sim)
    }
    return(sim)
}
