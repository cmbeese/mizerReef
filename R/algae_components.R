#' Algae Biomass
#'
#' The algae resource pool represents a combination of algal turf mats,
#' macroalgae, and the epilithic algal matrix (not including detritus).
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

# #' Algae dynamics with carrying capacity
# #'
# #' Calculates the algae biomass at the next time step from the current
# #' algae biomass
# #'
# #' The time evolution of the algae biomass \eqn{B} is described by
# #'
# #'  \deqn{ \frac{dB}{dt} = P_A\left( 1 - \frac{B}{\Kappa} \right) - c_A \, B_A }{
# #'                 dB/dt = P_A * (1 - B_A/\Kappa) - c_a * B_A}
# #'
# #' where \eqn{\Kappa} is the system's carrying capacity for algae in grams, 
# #' \eqn{c_A} is the mass-specific rate of consumption calculated with 
# #' `algae_consumption()` and \eqn{P_A} is the rate at which algae 
# #' grows, calculated with `getAlgaeProduction()`.
# #'
# #' The dynamical equation is solved analytically to
# #'       
# #'   \deqn{ B_A(t + dt) = B_A(t) \cdot e^{-\frac{dt}{\Kappa}(P_A+\Kappa \, C_A)} 
# #'                        - \frac{\Kappa \, P_A}{P_A+\Kappa \, c_A} 
# #'                        e^{-\frac{dt}{\Kappa}(P_A+\Kappa \, C_A)} }{
# #'         B_A(t + dt) = B_A(t) exp^(-dt/K * (P_A+ k*C_A))
# #'                        (-k*P_A) / (P_A + K*c_A) * e^(-dt/k * (P_A + k*C_A) }
# #'                        
# #' @param params A [MizerParams] object
# #' @param n A matrix of current species abundances (species x size)
# #' @param n_other Other dynamic components.
# #' @param rates A list of rates as returned by [getRates()]
# #' @param dt Time step size
# #' @param ... Unused
# #'
# #' @return A single number giving the algae biomass at next time step
# #' @seealso [detritus_dynamics()], [algae_consumption()], 
# #'          [getAlgaeConsumption()], [getAlgaeProduction()]
# #' @concept algae
# #' @export
# algae_dynamics2 <- function(params, n, n_other, rates, dt, ...) {
#     
#     consumption <- algae_consumption(params, n, rates)
#     production <- sum(getAlgaeProduction(params))
#     k <- params@other_params$algae$carry
#     
#     # If consumption is non-zero, return analytic solution
#     if (consumption) {
#         et <- exp(-dt/k * (production + k * consumption))
#         frac <- (k*production) / (production + k * consumption)
#         return(n_other$algae * et - frac * et)
#     }
#     return(n_other$algae + production * dt)
# }


#' Algae dynamics
#'
#' Calculates the algae biomass at the next time step from the current
#' algae biomass
#'
#' The time evolution of the algae biomass \eqn{B} is described by
#'
#' \deqn{dB/dt = \tt{production} - \tt{consumption} \cdot B}{
#'       dB/dt = production - consumption * B}
#'
#' where  `consumption` is the mass-specific rate of consumption calculated
#' with `algae_consumption()` and `production` is the rate at which algae 
#' grows, calculated with `getAlgaeProduction()`.
#'
#' The dynamical equation is solved analytically to
#'
#' \deqn{B(t+dt) = B(t) e^{(-\tt{consumption} \cdot dt)}
#'              +\frac{\tt{production}}{\tt{consumption}}
#'              (1- e^{(-\tt{consumption} \cdot dt)}).}{
#'       B(t+dt) = B(t) exp(-consumption * dt) 
#'               + production/consumption * (1 - exp(-consumption * dt)).}
#'
#' This avoids the stability problems that would arise if we used the Euler
#' method to solve the equation numerically.
#'
#' @param params A [MizerParams] object
#' @param n A matrix of current species abundances (species x size)
#' @param n_other Other dynamic components.
#' @param rates A list of rates as returned by [getRates()]
#' @param dt Time step size
#' @param ... Unused
#'
#' @return A single number giving the algae biomass at next time step
#' @seealso [detritus_dynamics()], [algae_consumption()], 
#'          [getAlgaeConsumption()], [getAlgaeProduction()]
#' @concept algae
#' @export
algae_dynamics <- function(params, n, n_other, rates, dt, ...) {

    consumption <- algae_consumption(params, n, rates)
    production <- sum(getAlgaeProduction(params))

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
#'  The rate at which herbivorous consumer groups encounter algal 
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

    # With feeding level
    feeding_level <- rates$feeding_level
    feeding_level[is.na(feeding_level)] <- 0
    sum((params@other_params$algae$rho * n * (1 - feeding_level))
        %*% params@dw)

    # # Without feeding level
    # sum((params@other_params$algae$rho * n) %*% params@dw)
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
    feeding_level[is.na(feeding_level)] <- 0
    consumption <- (params@other_params$algae$rho * params@initial_n *
                        (1 - feeding_level)) %*% params@dw

    # Without feeding levels
    # consumption <- (params@other_params$algae$rho * params@initial_n) %*% params@dw

    # Fix names
    names(consumption) <- params@species_params$species
    # Convert from mass specific rate to total rates
    consumption <- consumption * params@initial_n_other$algae

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
#' biomass. The rate is set so that production and consumption are equal for
#' chosen steady state abundances.
#'
#' @param params MizerParams
#'
#' @return The annual growth rate of algae per square meter
#' @seealso [getAlgaeConsumption()], [algae_dynamics()], [getAlgaeProduction()]
#' @concept algae
#' @export
getAlgaeProduction <- function(params) {
    return(params@other_params$algae$growth)
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


#' Expected algae lifetime
#'
#' The expected algae lifetime is defined as the inverse of the
#' mass-specific algae consumption rate.
#'
#' @param params A MizerParams object
#' @return The number giving the expected lifetime in years.
#' @concept algae
#' @export
algae_lifetime <- function(params) {
    1 / algae_consumption(params,
                          n = params@initial_n,
                          rates = getRates(params))
}

#' @rdname algae_lifetime
#' @param params A MizerParams object
#' @param value A number with the new value for the expected lifetime in years
#'
#' Assigning a new value to the algae lifetime rescales the algae
#' abundance while keeping the total consumption of algae the same (by
#' adjusting the interaction strength of species with algae).
#' 
#' @concept algae
#' @export
`algae_lifetime<-` <- function(params, value) {
    rescale_algae(params, value / algae_lifetime(params))
}

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
    params@other_params[["algae"]]$rho <-
        params@other_params[["algae"]]$rho / factor
    params
}
