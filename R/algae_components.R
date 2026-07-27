#' algae Biomass
#'
#' The algae resource pool represents a combination of algae turf mats,
#' macroalgae, and the epilithic algae matrix (not including detritus).
#' It is not size structured to reflect the fact that herbivorous fish
#' feed on algae regardless of their body size.
#'
#' @param params MizerParams
#' @return The algae biomass in grams
#' @concept algae
#' @export
algae_biomass <- function(params) {
    params@initial_n_other$algae
}

#' Algae dynamics with carrying capacity
#'
#' Calculates the algae biomass at the next time step from the current
#' algae biomass
#'
#' The time evolution of the algae biomass \eqn{B_A(t)} is described by
#'
#'  \deqn{ \frac{dB_A}{dt} = P_A\left( 1 - 
#'                          \frac{B_A}{K_A} \right) - c_A \, B_A }{
#'                 dB_A/dt = P_A * (1 - B_A/ K_A) - c_A * B_A}
#'
#' where \eqn{K_A} is the system's carrying capacity for algae in grams/ year,
#' \eqn{c_A} is the mass-specific rate of consumption calculated with
#' `algae_consumption()` and \eqn{P_A} is the rate at which algae
#' grows, calculated with `getAlgaeProduction()`.
#'
#' The dynamical equation is solved analytically to
#'
#'   \deqn{B_A(t + dt) = B_A(t) \cdot e^{-\frac{dt}{K_A}(P_A+ K_A \, c_A)}
#'                        + \frac{K_A \, P_A}{P_A+ K_A \, c_A} \left(1-
#'                        e^{-\frac{dt}{K_A}(P_A+ K_A \, c_A)}\right) }{
#'         B_A(t + dt) = B_A(t) exp^(-dt/K_A * (P_A+ K_A*c_A)) +
#'                        (K_A*P_A) / (P_A + K_A*c_A) *(1 - e^(-dt/K_A * 
#'                        (P_A + K_A*c_A)) }
#'
#' @param params A [MizerParams] object
#' @param n A matrix of current species abundances (species x size)
#' @param n_other Other dynamic components.
#' @param rates A list of rates as returned by [getRates()]
#' @param dt Time step size
#' @param t The current time, used to determine the post-bleaching algae
#'          growth and capacity boosts, if any (see [getAlgaeBoost()]).
#'          Defaults to `0`, which is always before any bleaching event.
#' @param ... Unused
#'
#' @return A single number giving the algae biomass at next time step
#' @seealso [detritus_dynamics()], [algae_consumption()],
#'          [getAlgaeConsumption()], [getAlgaeProduction()], [getAlgaeBoost()]
#' @concept algae
#' @export
algae_dynamics_cc <- function(params, n, n_other, rates, dt, t = 0, ...) {

    consumption <- algae_consumption(params, n, rates)
    production <- sum(getAlgaeProduction(params, t))
    capacity_boost <- getAlgaeBoost(
        params, t,
        params@other_params$algae_params$algae_capacity_boost
    )
    ka <- params@other_params$algae_params$algae_capacity * capacity_boost

    if(is.nan(consumption)){ 
        warning("The algae consumption function is producing NaNs.")
        }

    # If consumption is non-zero, return analytic solution
    if (consumption) {
        et <- exp(-dt/ka * (production + ka * consumption))
        frac <- (ka*production) / (production + ka * consumption)
        fracet <- frac *(1- et)
        return(n_other$algae * et + fracet)
    }
    et <- exp(-dt/ka * (production))
    return(n_other$algae * et)
}


#' Algae dynamics
#'
#' Calculates the algal biomass at the next time step from the current
#' algae biomass
#'
#' The time evolution of the algal biomass \eqn{B} is described by
#'
#' \deqn{dB_A/dt = P_A - c_A \cdot B_A}{
#'       dB_A/dt = P_A - c_A * B_A}
#'
#' where  \eqn{c_A} is the mass-specific rate of consumption calculated
#' with `algae_consumption()` and \eqn{P_A} is the rate at which algae 
#' grows, calculated with `getAlgaeProduction()`.
#'
#' The dynamical equation is solved analytically to
#'
#' \deqn{B_A(t+dt) = B_A(t) e^{(- c_A \cdot dt)}
#'              +\frac{P_A}{c_A}
#'              (1- e^{(-c_A \cdot dt)}).}{
#'       B_A(t+dt) = B(t) exp(-c_A * dt) 
#'               + p_A/c_A * (1 - exp(-c_A * dt)).}
#'
#' This avoids the stability problems that would arise if we used the Euler
#' method to solve the equation numerically.
#'
#' @param params A [MizerParams] object
#' @param n A matrix of current species abundances (species x size)
#' @param n_other Other dynamic components.
#' @param rates A list of rates as returned by [getRates()]
#' @param dt Time step size
#' @param t The current time, used to determine the post-bleaching algae
#'          growth boost, if any (see [getAlgaeBoost()]). Defaults to `0`,
#'          which is always before any bleaching event.
#' @param ... Unused
#'
#' @return A single number giving the algae biomass at next time step
#' @seealso [detritus_dynamics()], [algae_consumption()],
#'          [getAlgaeConsumption()], [getAlgaeProduction()], [getAlgaeBoost()]
#' @concept algae
#' @export
algae_dynamics <- function(params, n, n_other, rates, dt, t = 0, ...) {

    consumption <- algae_consumption(params, n, rates)
    production  <- sum(getAlgaeProduction(params, t))
    
    if(is.nan(consumption)){ 
        warning("The algae consumption function is producing NaNs.")
    }
    

    # If consumption is non-zero, return analytic solution
    if (consumption) {
        et <- exp(-consumption * dt)
        return(n_other$algae * et + production / consumption  * (1 - et))
    }
    return(n_other$algae + production * dt)
}

#' Mass-specific algae consumption rate
#'
#' This mass-specific consumption rate is used in `algae_dynamics()` to
#' calculate the algae biomass at the next time step. To get the
#' non-mass-specific consumption rate, use `getAlgaeConsumption()`.
#' 
#' The rho parameter for herbivorous fish groups is stored in
#' `other_params(params)$algae$rho`
#'
#' @section Algae consumption:
#' 
#'  The rate at which herbivorous consumer groups encounter algae 
#'  biomass \eqn{E_{i.A}(w)} is controlled by the parameter 
#'  \eqn{\rho_{A.i}}. It scales with the size of the consumer raised to 
#'  an allometric exponent \eqn{m_{alg}} which is taken from empirical data.
#'  
#'  \deqn{E_{i.A}(w)=\rho_{i.A}\, w^{m_{alg}}\,B_A}{
#'        E_{i.A}(w)= rho_{i.A}\, w^{m_{alg}}\,B_A}
#'          
#'  The mass specific consumption rate then accounts for the preference of 
#'  functional group $i$ for algae, \eqn{\theta_{i.A}}. This gives the 
#'  mass-specific algae consumption rate:
#'  
#'  \deqn{c_A = \sum_i\int\rho_{i.A}\, w^{m_{alg}} N_i(w)\theta_{i.A}\,dw}{
#'        c_A = \sum_i\int rho_{i.A}\, w^{m_{alg}} N_i(w) theta_{i.A}\,dw}
#'
#' @param params MizerParams
#' @param n A matrix of current species abundances (species x size)
#' @param rates A list of rates as returned by [getRates()]
#'
#' @return The mass-specific consumption rate of algae in grams per year.
#' @concept algae
#' @export
algae_consumption <- function(params,
                              n = params@initial_n,
                              rates = getRates(params)) {

    sum((params@other_params$algae_params$rho * n) %*% params@dw)
}

#' Get algae consumption rates
#'
#' This function returns a named vector with one component for each species
#' giving the rate in grams/year at which that species consumes algae
#' 
#' @inheritSection algae_consumption Algae consumption
#'
#' @param params MizerParams
#' @return A named vector with the consumption rates from herbivores
#' @seealso [getAlgaeProduction()], [algae_dynamics()], [getAlgaeConsumption()]
#' @concept algae
#' @export
getAlgaeConsumption <- function(params) {

    # With feeding level
    feeding_level <- getFeedingLevel(params)
    consumption <- (params@other_params$algae_params$rho * params@initial_n *
                        (1 - feeding_level)) %*% params@dw

    # Convert from mass specific rate to total rates
    consumption <- consumption * params@initial_n_other$algae

    # Fix names
    names(consumption) <- params@species_params$species

    return(consumption)
}

#' Plot algae consumption rates
#'
#' @param params MizerParams
#' @return A pie chart.
#' @concept algae
#' @export
plotAlgaeConsumption <- function(params) {
    consumption <- getAlgaeConsumption(params)
    total <- sum(consumption)
    consumption <- consumption[consumption > total/100]
    df <- data.frame(Consumer = names(consumption),
                     Rate = consumption)
    ggplot(df, aes(x = "", y = Rate, fill = Consumer)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        labs(title = "Algae consumption rate [g/year]",
             x = "", y = "")
}

#' Algae production rate
#'
#' This is the rate in grams/year/m^-2 at which the system produces algae
#' biomass. The baseline rate is set so that production and consumption are
#' equal for chosen steady state abundances (see [tuneUR()]/[tuneUR_cc()]).
#' If post-bleaching algae boosting has been configured with
#' [setDegradation()], the baseline rate is scaled by the boost factor
#' appropriate for time `t` (see [getAlgaeBoost()]).
#'
#' @param params MizerParams
#' @param t The current time, used to determine the post-bleaching algae
#'          growth boost, if any (see [getAlgaeBoost()]). Defaults to `0`,
#'          which is always before any bleaching event, so no boost applies.
#'
#' @return The annual growth rate of algae per square meter
#' @seealso [getAlgaeConsumption()], [algae_dynamics()], [getAlgaeBoost()]
#' @concept algae
#' @export
getAlgaeProduction <- function(params, t = 0) {
    boost <- getAlgaeBoost(
        params, t,
        params@other_params$algae_params$algae_growth_boost
    )
    params@other_params$algae_params$algae_growth * boost
}

#' Post-bleaching boost multiplier for algae growth or capacity
#'
#' Computes the multiplier that [getAlgaeProduction()] and
#' [algae_dynamics_cc()] apply to the baseline algae growth rate or carrying
#' capacity to represent a temporary post-bleaching boost in algal
#' productivity, as configured by the `algae_boost`, `algae_growth_boost` and
#' `algae_capacity_boost` arguments to [setDegradation()].
#'
#' This is a pure function of `t`: it returns the cumulative product of
#' `boost_vector` up to the element for the current year post-bleaching,
#' applied fresh to the untouched baseline value every time it is called.
#' This mirrors how [reefDegrade()] recomputes refuge density fresh from the
#' baseline and `deg_scale` at every timestep, rather than by mutating
#' `params`: because a single `mizer::project()` call treats `params` as
#' read-only across the whole simulation (only `n`/`n_pp`/`n_other` are
#' threaded forward between timesteps), a rate function cannot rely on a
#' mutation to `params` persisting into the next timestep.
#'
#' @param params A MizerParams object
#' @param t The current time
#' @param boost_vector The boost vector to use -- `algae_growth_boost` or
#'   `algae_capacity_boost` from `params@other_params$algae_params`.
#'   Element 1 gives the multiplier for the bleaching year itself, element 2
#'   for one year post-bleaching, and so on. If `t` is later than the last
#'   element of `boost_vector`, the cumulative product stays fixed at the
#'   product of the whole vector (i.e. the boost plateaus rather than
#'   continuing to compound).
#'
#' @return A single number: `1` (no boost) if `algae_boost` is not `TRUE`,
#'   `boost_vector` is empty, or `t` is before the bleaching year; otherwise
#'   the cumulative product of `boost_vector` up to the current post-bleach
#'   year.
#' @seealso [getAlgaeProduction()], [algae_dynamics_cc()], [setDegradation()],
#'   [reefDegrade()]
#' @concept algae
#' @export
getAlgaeBoost <- function(params, t, boost_vector) {
    if (!isTRUE(params@other_params$algae_params$algae_boost) ||
        is.null(boost_vector) || length(boost_vector) == 0) {
        return(1)
    }
    t_bleach <- params@other_params$refuge_params$t_bleach
    if (is.null(t_bleach) || t < t_bleach) {
        return(1)
    }
    years_post <- t - t_bleach
    idx <- min(years_post + 1, length(boost_vector))
    prod(boost_vector[seq_len(idx)])
}


# Probably not a useful function at the moment since algae production 
# is constant (for now)
# #' Plot algae production rates
# #'
# #' @param params MizerParams
# #' @return A pie chart.
# #' @export
# plotalgaeProduction <- function(params) {
#   production <- getalgaeProduction(params)
#   df <- data.frame(Producer = names(production),
#                    Rate = production)
#   ggplot(df, aes(x = "", y = Rate, fill = Producer)) +
#     geom_bar(stat = "identity", width = 1) +
#     coord_polar("y", start = 0) +
#     labs(title = "algae production rate [g/year]",
#          x = "", y = "")
# }

# Not used
# #' Expected algae lifetime
# #'
# #' The expected algae lifetime is defined as the inverse of the
# #' mass-specific algae consumption rate.
# #'
# #' @param params A MizerParams object
# #' @return The number giving the expected lifetime in years.
# #' @concept algae
# #' @export
# algae_lifetime <- function(params) {
#     1 / algae_consumption(params,
#                           n = params@initial_n,
#                           rates = getRates(params))
# }
# 
# #' @rdname algae_lifetime
# #' @param params A MizerParams object
# #' @param value A number with the new value for the expected lifetime in years
# #'
# #' Assigning a new value to the algae lifetime rescales the algae
# #' abundance while keeping the total consumption of algae the same (by
# #' adjusting the interaction strength of species with algae).
# #' 
# #' @concept algae
# #' @export
# `algae_lifetime<-` <- function(params, value) {
#     rescale_algae(params, value / algae_lifetime(params))
# }

#' Rescale algae biomass without changing anything else
#'
#' This multiplies the algae biomass by a factor and divides the
#' interaction between all species and algae by the same
#' factor, so as to keep the total consumption of algae unchanged.
#'
#' @param params A MizerParams object
#' @param factor A number
#' @return An updated MizerParams object
#' @concept algae
#' @export
rescale_algae <- function(params, factor) {
    params@initial_n_other[["algae"]] <-
        params@initial_n_other[["algae"]] * factor
    params@species_params$rho_algae <-
        params@species_params$rho_algae / factor
    params@other_params$algae_params$rho <-
        params@other_params$algae_params$rho / factor
    params
}
