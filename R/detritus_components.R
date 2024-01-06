#' Detritus Biomass
#'
#' The detrital resource pool represents any consumed wastes, including
#' decomposing dead organisms and feces. Similar to algae, the detritus
#' resource is not size structured because herbivores of any size feed on
#' detritus on coral reefs.
#'
#' @param params MizerParams
#' @return The detritus biomass in grams
#' @concept detritus
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
#' \deqn{dB/dt = \tt{production} - \tt{consumption} \cdot B + \tt{external}}{
#'       dB/dt = production - consumption * B + external}
#'
#' where `consumption` is the mass-specific rate of consumption, calculated 
#' with `detritus_consumption()`and `production` is the rate at which the 
#' rest of the system produces detritus biomass, calculate with 
#' `getDetritusProduction()`.
#' 
#' The dynamical equation is solved analytically to
#' \deqn{B(t+dt) = B(t)\exp(-\tt{consumption} \cdot dt)
#'               + \frac{\tt{production}}{\tt{consumption}}
#'                 (1-\exp(-\tt{consumption} \cdot dt)).}{
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
#' @return A vector giving the detritus spectrum at the next time step.
#' @seealso [algae_dynamics()], [detritus_consumption()], 
#'          [getDetritusConsumption()], [getDetritusProduction()]
#' @concept detritus
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
#' The rho parameter for each functional group is stored in
#' `other_params(params)$detritus$rho`
#' 
#' @inheritSection getDetritusConsumption Detritus consumption
#'
#' @param params MizerParams
#' @param n A matrix of current species abundances (species x size)
#' @param rates A list of rates as returned by [getRates()]
#'
#' @return The mass-specific consumption rate of detritus in grams per year.
#' @concept detritus
#' @export
detritus_consumption <- function(params,
                                 n = params@initial_n,
                                 rates = getRates(params)) {
    # With feeding level
    feeding_level <- rates$feeding_level
    feeding_level[is.na(feeding_level)] <- 0
    sum((params@other_params$detritus$rho * n * (1 - feeding_level))
        %*% params@dw)

    # Without feeding level
    #sum((params@other_params$detritus$rho * n ) %*% params@dw)
}


#' Get detritus consumption rates
#'
#' This function returns a named vector with one component for each species
#' giving the rate in grams/year at which that species consumes detritus
#' 
#' @section Detritus consumption:
#' 
#'  The rate at which detritivorous consumer groups encounter detrital 
#'  biomass \eqn{E_{i.D}(w)} is controlled by the parameter 
#'  \eqn{\rho_{D.i}}. It scales with the size of the consumer raised to 
#'  an allometric exponent \eqn{m_{det}} which is taken to be the same as 
#'  the scaling exponent of the maximum intake rate for fish consumers.
#'  
#'  \deqn{E_{i.D}(w)=\rho_{i.D}\, w^{m_{det}}\,B_D }{
#'        E_{i.D}(w)=\rho_{i.D}\, w^{m_{det}}\,B_D}
#'          
#'  The mass specific consumption rate then accounts for the preference of 
#'  functional group $i$ for detritus, \eqn{\theta_{i.D}} and the feeding 
#'  level \eqn{f_i(w)}. This gives the mass-specific detritus consumption
#'  rate:
#'  
#'  \deqn{c_D = \sum_i\int\rho_{i.D}\, w^{m_{det}} 
#'              N_i(w) \left(1-f_i(w)\right) \theta_{i.D}\,dw}{
#'        c_D = \sum_i\int\rho_{i.D}\, w^{m_{det}} 
#'              N_i(w) \left(1-f_i(w)\right) \theta_{i.D}\,dw}
#'              
#' 
#' @param params MizerParams
#' @return A named vector with the consumption rates from herbivores
#' @seealso [getAlgaeProduction()], [algae_dynamics()], [getDetritusConsumption()]
#' 
#' @concept detritus
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
#' @concept detritus
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
#' \enumerate{
#'      \item consumed biomass not assimilated by predators ("feces"),
#'      \item decomposing dead organisms ("decomp"),
#'      \item the pelagic zone ("external").
#'  }
#'  
#' This function returns a vector with the individual contributions for each
#' source. These can be summed with `sum()` to get the total detritus
#' production rate.
#' 
#' @section Detritus production:
#'    
#'  The rate \eqn{p_D} at which detritus biomass is produced by the 
#'  ecosystem has contributions from three sources:
#'  
#'  \deqn{p_D = p_{D.f} + p_{D.d} + p_{D.ext}}{
#'        p_D = p_{D.f} + p_{D.d} + p_{D.ext}}
#'        
#'  \eqn{p_{D.f}} comes from the biomass that is consumed but not 
#'  assimilated and is given by:
#'  
#'  \deqn{p_{D.f} = \sum_i(1-\alpha_i)\int E_i(w)\,dw}{
#'        p_{D.f} = \sum_i(1-\alpha_i)\int E_i(w)\,dw}
#'        
#'  \eqn{p_{D.d}} comes from the biomass of fish that die as a result of 
#'  external mortality. External mortality includes local deaths that lead 
#'  to detritus but also deaths due to predation by species that are not
#'  explicitly modelled, for example transient predators, mammals, or sea 
#'  birds. Thus, only a proportion `prop_decomp` of this material 
#'  decomposes to detritus. The detritus production from decomposing 
#'  dead organisms is given by:
#'  
#'  \deqn{p_{D.d} = \sum_i\int\mu_{seni.i}(w)N_i(w)w\,dw + 
#'                  \mathtt{prop\_decomp}\,
#'                  \sum_i\int\mu_{nat.i}(w)N_i(w)w\,dw}{
#'      p_{D.d} = \sum_i\int\mu_{seni.i}(w)N_i(w)w\,dw + 
#'                  \mathtt{prop\_decomp}\,
#'                  \sum_i\int\mu_{nat.i}(w)N_i(w)w\,dw}
#'                  
#'  \eqn{p_{D.ext}} is the rate at which detritus enters the system from
#'  unmodelled or external sources. For coral reefs, this includes detritus
#'  produced by sponges and coral mucous as well as waste material that 
#'  sinks in from the pelagic zone. This rate is a model parameter 
#'  independent of any other model component. It is set so that production 
#'  and consumption are equal for the chosen steady state abundances.
#'
#' @param params MizerParams
#' @param n A matrix of current species abundances (species x size)
#' @param rates A list of rates as returned by [getRates()]
#'
#' @return A vector with named entries "feces", "decomp", and "external",
#' giving the rates at which detritus biomass is produced by each of these
#' sources in grams per year.
#' @concept detritus
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
#' @concept detritus
#' @export
plotDetritusProduction <- function(params) {
    production <- getDetritusProduction(params)
    df <- data.frame(Source = names(production),
                     Rate = production)
    ggplot(df, aes(x = "", y = Rate, fill = Source)) +
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
#' @concept detritus
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
#' @concept detritus
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
#' @concept detritus
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
