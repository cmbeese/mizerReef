#' Contribution of unstructured components to the encounter rate
#'
#' The encounter rate \eqn{E_i(w)} for an unstructured resource like algae or
#' detritus is proportional to the total biomass \eqn{B} with a coefficient
#' \eqn{\rho_i(w)} that depends on the functional group \eqn{i} and the size of
#' the consumer:
#'
#' \deqn{E_i(w) = \rho_i(w) B.}
#'
#' The coefficient \eqn{\rho_i(w)} is stored as a matrix (species x size) in
#' the `rho` parameter of the component. It has units 1/year.
#'
#' @param params MizerParams
#' @param n_other Biomasses of unstructured components
#' @param component Name of component whose contribution is requested
#' @param ... Unused
#'
#' @return Array (species x size) with the encounter rate in g/year.
#' @export
#' @concept Uresources
encounter_contribution <- function(params, n_other, component, ...) {
    
    params@other_params[[component]]$rho * n_other[[component]]
}

#' Rescale algae and detritus biomass without changing anything else
#'
#' This multiplies the algae & detritus biomass by a factor and divides the
#' interaction between all groups and unstructured resource by the same
#' factor, so as to keep the total consumption of these resources unchanged.
#' 
#' @param params A MizerParams object
#' @param algae_factor A number to scale algae by
#' @param detritus_factor A number to scale detritus by
#' @return An updated MizerParams object
#' @concept Uresources
#' @export
rescaleComponents <- function(params, algae_factor = 1, detritus_factor = 1) {
    rescale_algae(rescale_detritus(params, detritus_factor),
                  algae_factor)
}

#' Tune unstructured resources (algae and detritus) to steady state
#'
#' This first sets the rate of degradation of algae so that for the given
#' abundances, the algae is at steady state. It then sets the rate at which
#' detritus flows in from external sources (e.g. the pelagic zone) so that for
#' the given abundances the detritus is at steady state.
#'
#' @param params A MizerParams object
#' @param ... unused
#' @return An updated MizerParams object
#' @concept Uresources
#' 
#' @export
#' 
tuneUR <- function(params,...) {

    # algae
    ain  <- getAlgaeProduction(params) / params@initial_n_other$algae
    aout <- algae_consumption(params)
    if (ain < aout) {
        warning("The value for algae growth provided does not produce enough
                to support this abundance of herbivores. I will increase algae
                growth to meet the consumption rate.")
    }
    params@other_params$algae$algae_growth <- aout

    # detritus
    params@other_params$detritus$external <- 0
    din <- sum(getDetritusProduction(params))
    dout <- detritus_consumption(params)
    if (din < dout) {
        warning("Detrital production is not high enough to support this 
                abundance of detritivores. I will increase external
                input to meet the consumption rate.")
    }
    params@other_params$detritus$external <- dout - din

    params
}

#' Scale Model Parameters
#'
#' This function scales various model parameters by a given factor.
#'
#' @param params a mizer model object
#' @param factor a numeric value by which to scale the model
#'
#' @return a mizer model object with scaled parameters
#' @concept Uresources
#' @export
reef_scaleModel <- function(params, factor) {

    # Algae
    params@other_params[["algae"]]$rho <-
        params@other_params[["algae"]]$rho / factor
    params@species_params$rho_algae <-
        params@species_params$rho_algae / factor

    # Detritus
    params@other_params[["detritus"]]$rho <-
        params@other_params[["detritus"]]$rho / factor
    params@species_params$rho_detritus <-
        params@species_params$rho_detritus / factor
    params@other_params$detritus$external <-
        params@other_params$detritus$external * factor

    # now comes the code of mizer's standard scaleModel()
    params <- validParams(params)
    assert_that(is.number(factor), factor > 0)
    params@cc_pp <- params@cc_pp * factor
    params@resource_params$kappa <- params@resource_params$kappa * 
        factor
    if ("r_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$r_max
        params@species_params$r_max <- NULL
        message("The 'r_max' column has been renamed to 'R_max'.")
    }
    if ("R_max" %in% names(params@species_params)) {
        params@species_params$R_max <- params@species_params$R_max * 
            factor
    }
    params@search_vol <- params@search_vol/factor
    if ("gamma" %in% names(params@species_params)) {
        params@species_params$gamma <- params@species_params$gamma/factor
    }
    initial_n_other <- params@initial_n_other
    for (res in names(initial_n_other)) {
        initial_n_other[[res]] <- initial_n_other[[res]] * factor
    }
    initialN(params) <- params@initial_n * factor
    initialNResource(params) <- params@initial_n_pp * factor
    initialNOther(params) <- initial_n_other
    params@sc <- params@sc * factor
    return(params)
}

#' Hold resource dynamics constant
#' 
#' @param params MizerParams object
#' @param n_other Biomasses of unstructured components
#' @param component Name of component to view dynamics for
#' @param ... Unused
#' @concept Uresources
#' @export
constant_dynamics <- function(params, n_other, component, ...) {
    n_other[[component]]
}
