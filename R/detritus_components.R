#' Detritus Biomass
#'
#' The detrital resource pool represents any consumed wastes, including
#' decomposing dead organisms and feces. Similar to algae, the detritus
#' resource is not size structured because herbivores of any size feed on
#' detritus on coral reefs.
#'
#' @param params MizerParams
#' @return The detritus biomass in grams
#' @concept unstructured resources
#' @export
detritus_biomass <- function(params) {
    params@initial_n_other$detritus
}


#' Detritus dynamics
#'
#' Calculates the detritus biomass at the next time step based on the
#' current detritus biomass.
#'
#' The time evolution of the detritus biomass \eqn{B} is described by
#'
#' \deqn{dB/dt = \tt{production} - \tt{consumption} * B + \tt{external}}{dB/dt = production - consumption * B + external}
#'
#' where
#'  * `consumption` is the mass-specific rate of consumption, calculated with
#'     `detritus_consumption()`
#'  * `production` is the rate at which the rest of the system produces
#'    detritus biomass, calculate with `getDetritusProduction()`
#'
#' The dynamical equation is solved analytically to
#' \deqn{B(t+dt) = B(t)\exp(-\tt{consumption} \cdot dt)
#'   +\frac{\tt{production}}{\tt{consumption}}
#'   (1-\exp(-\tt{consumption} \cdot dt)).}{B(t+dt)
#'   = B(t) exp(-consumption * dt)
#'   +production/consumption * (1 - exp(-consumption * dt)).}
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
#' @return A vector giving the detritus spectrum at the next time step.
#' @seealso [algae_dynamics()]
#' @concept unstructured resources
#' @export
detritus_dynamics <- function(params, n, n_other, rates, dt, ...) {

    consumption <- detritus_consumption(params, n, rates)
    production <- sum(getDetritusProduction(params, n, rates))

    if (consumption) {
        et <- exp(-consumption * dt)
        return(n_other$detritus * et + production / consumption  * (1 - et))
    }
    return(n_other$detritus + production * dt)
}

#' Mass-specific detritus consumption rate
#'
#' This mass-specific consumption rate is used in `detritus_dynamics()` to
#' calculate the detritus biomass at the next time step. To get the
#' non-mass-specific consumption rate, use `getDetritusConsumption()`.
#'
#' The consumption rate by detritivorous fish is determined by
#' `other_params(params)$detritus$rho`
#'
#' @param params MizerParams
#' @param n A matrix of current species abundances (species x size)
#' @param rates A list of rates as returned by [getRates()]
#'
#' @return The mass-specific consumption rate of detritus in grams per year.
#' @concept unstructured resources
#' @export
detritus_consumption <- function(params,
                                 n = params@initial_n,
                                 rates = getRates(params)) {
    # With feeding level
    sum((params@other_params$detritus$rho * n * (1 - rates$feeding_level))
        %*% params@dw)

    # Without feeding level
    #sum((params@other_params$detritus$rho * n ) %*% params@dw)
}


#' Get detritus consumption rates
#'
#' This function returns a named vector with one component for each species
#' giving the rate in grams/year at which that species consumes detritus
#' @param params MizerParams
#' @return A named vector with the consumption rates from herbivores
#' @seealso [getalgaeProduction()], [algae_dynamics()], [getDetritusConsumption()]
#' 
#' @concept unstructured resources
#' @export
getDetritusConsumption <- function(params) {

    # With feeding level
    feeding_level <- getFeedingLevel(params)
    consumption <- (params@other_params$detritus$rho * params@initial_n *
                        (1 - feeding_level)) %*% params@dw

    # Without feeding level
    # consumption <- (params@other_params$detritus$rho * params@initial_n) %*% params@dw

    names(consumption) <- params@species_params$species
    # Convert from mass specific rate to total rates
    consumption <- consumption * params@initial_n_other$detritus

    return(consumption)
}


#' Plot detritus consumption rates by species
#'
#' @param params MizerParams
#' @return A pie chart.
#' @concept unstructured resources
#' @export
plotDetritusConsumption <- function(params) {
    consumption <- getDetritusConsumption(params)
    total <- sum(consumption)
    consumption <- consumption[consumption > total/100]
    df <- data.frame(Consumer = names(consumption),
                     Rate = consumption)
    ggplot(df, aes(x = "", y = Rate, fill = Consumer)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        labs(title = "Detritus consumption rate [g/year]",
             x = "", y = "")
}


#' Detritus production rate
#'
#' Gives a named vector with the rates at which different components of the
#' ecosystem produce detritus:
#'
#' 1. consumed biomass not assimilated by predators ("feces"),
#' 2. decomposing dead organisms ("decomp"),
#' 3. the pelagic zone ("external").
#'
#' mizerReef models include two sources of external mortality to describe all
#' deaths by natural causes that are not due to predation by the modelled
#' species. External mortality includes deaths that lead to detritus, but also
#' deaths due to predation by species that are not explicitly modelled, for
#' example transient predators, mammals, or sea birds. Thus, only a proportion
#' of this material will stay in the system to decompose. The parameter
#' `prop_decomp` describes the proportion of external mortality that stays in
#' the system and decomposes to detritus.
#'
#' External detritus is waste material that sinks in from the pelagic zone.
#' This rate is a model parameter independent of any other model component. It
#' is set so that production and consumption are equal for chosen steady state
#' abundances.
#'
#' The function returns a vector with the individual contributions for each
#' source. These can be summed with `sum()` to get the total detritus
#' production rate.
#'
#' @param params MizerParams
#' @param n A matrix of current species abundances (species x size)
#' @param n_other Other dynamic components. Only `n_other$algae` is used.
#' @param rates A list of rates as returned by [getRates()]
#'
#' @return A vector with named entries "feces", "decomp", and "external",
#' giving the rates at which detritus biomass is produced by each of these
#' sources in grams per year.
#' @concept unstructured resources
#' @export
getDetritusProduction <- function(params, n = params@initial_n,
                                  rates = getRates(params)) {
    # Feces
    # With feeding level
    consumption <- sweep((1 - rates$feeding_level) * rates$encounter * n, 2,
                         params@dw, "*", check.margin = FALSE)

    # Without feeding level
    # consumption <- sweep(rates$encounter * n, 2,
    #                      params@dw, "*", check.margin = FALSE)

    feces <- sweep(consumption, 1, (1 - params@species_params$alpha), "*",
                   check.margin = FALSE)

    # Decomposition of dead organisms
    ex_mort <- sum((params@mu_b * n) %*% (params@w * params@dw))
    sen_mort <- getSenMort(params)
    sen_mort <- sum((sen_mort * n) %*% (params@w * params@dw))
    prop_decomp <- params@other_params$detritus$prop_decomp

    # Return vector
    c(feces    = sum(feces),
      decomp   = (prop_decomp*ex_mort) + sen_mort,
      external = params@other_params$detritus$d.external
    )
}


#' Plot detritus production rates from each source
#'
#' @param params MizerParams
#' @return A pie chart.
#' @concept plots
#' @export
plotDetritusProduction <- function(params) {
    production <- getDetritusProduction(params)
    df <- data.frame(Producer = names(production),
                     Rate = production)
    ggplot(df, aes(x = "", y = Rate, fill = Producer)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y", start = 0) +
        labs(title = "Detritus production rate [g/year]",
             x = "", y = "")
}


#' Expected detritus lifetime
#'
#' The expected detritus lifetime is defined as the inverse of the
#' mass-specific detritus consumption rate.
#'
#' @param params A MizerParams object
#' @return The number giving the expected lifetime in years.
#' @concept unstructured resources
#' @export
detritus_lifetime <- function(params) {
    1 / detritus_consumption(params,
                             n = params@initial_n,
                             rates = getRates(params))
}

#' @rdname detritus_lifetime
#' @param params A MizerParams object
#' @param value A number with the new value for the expected lifetime in years
#'
#' Assigning a new value to the detritus lifetime rescales the detritus
#' abundance while keeping the total consumption of detritus the same (by
#' adjusting the interaction strength of species with detritus).
#'
#' @concept unstructured resources
#' @export
`detritus_lifetime<-` <- function(params, value) {
    rescale_detritus(params, value / detritus_lifetime(params))
}

#' Rescale detritus biomass without changing anything else
#'
#' This multiplies the detritus biomass by a factor and divides the
#' interaction between all species and the detritus by the same
#' factor, so as to keep the total consumption of detritus unchanged.
#' It also divides the mass-specific rate of decomposition by the same
#' factor so that the total detritus decomposition rate stays the same.
#' @param params A MizerParams object
#' @param factor A number
#' 
#' @return An updated MizerParams object
#' 
#' @concept unstructured resources
#' @export
rescale_detritus <- function(params, factor) {
    params@initial_n_other[["detritus"]] <-
        params@initial_n_other[["detritus"]] * factor
    params@species_params$rho_detritus <-
        params@species_params$rho_detritus / factor
    params@other_params[["detritus"]]$rho <-
        params@other_params[["detritus"]]$rho / factor
    params
}
