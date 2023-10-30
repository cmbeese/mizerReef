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
#'
encounter_contribution <- function(params, n_other, component, ...) {
    params@other_params[[component]]$rho * n_other[[component]]
}


#' @export
rescaleComponents <- function(params, algae_factor = 1, detritus_factor = 1) {
    rescale_algae(rescale_detritus(params, detritus_factor),
                  algae_factor)
}

#' Tune algae and detritus to steady state
#'
#' This first sets the rate of degradation of algae so that for the given
#' abundances, the algae is at steady state. It then sets the rate at which
#' detritus flows in from external sources (e.g. the pelagic zone) so that for
#' the given abundances the detritus is at steady state.
#'
#' @param params A MizerParams object
#' @return An updated MizerParams object
#' @export
tune_algae_detritus <- function(params) {

    # algae
    ain <- sum(getAlgaeProduction(params)) / params@initial_n_other$algae
    aout <- algae_consumption(params)
    if (ain < aout) {
        stop("There is not enough algae production to support this abundance
             of herbivores.")
    }
    params@other_params$algae$algae_growth <- ain - aout

    # detritus
    params@other_params$detritus$external <- 0
    din <- sum(getDetritusProduction(params))
    dout <- detritus_consumption(params)
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
#'
#' @export
#'
#' @examples
#' params <- scaleModel(CBN_params, 0.5)
#'
scaleModel <- function(params, factor) {

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

    # Rest of mizer
    mizer::scaleModel(params, factor)
}

#' @export
constant_dynamics <- function(params, n_other, component, ...) {
    n_other[[component]]
}
