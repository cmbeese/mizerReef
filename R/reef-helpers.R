#' Remove some species from the model
#'
#' This extends mizer's [removeSpecies()] to also remove the relevant
#' row from the detritus and algae consumption arrays `rho`.
#' @param params A `mizerReef` object
#' @param species The species to be removed. A vector of species names, or a
#'   numeric vector of species indices, or a logical vector indicating for each
#'   species whether it is to be removed (TRUE) or not.
#' @param ... Unused
#' @return A `mizerReef` object with fewer species.
#' @method removeSpecies mizerReef
#' @export
#' @concept helper
removeSpecies.mizerReef <- function(params, species, ...) {
    keep <- !valid_species_arg(params, species, return.logical = TRUE)

    p <- NextMethod()

    # Remove algae rho values for species
    p@other_params$algae_params$rho <-
        p@other_params$algae_params$rho[keep, , drop = FALSE]

    # Remove detritus rho values for species
    p@other_params$detritus_params$rho <-
        p@other_params$detritus_params$rho[keep, , drop = FALSE]

    p
}

#' Check whether two objects are different
#'
#' Check whether two objects are numerically different, ignoring all attributes.
#'
#' We use this helper function in particular to see if a new value for a slot
#' in MizerParams is different from the existing value in order to give the
#' appropriate messages.
#'
#' @param a First object
#' @param b Second object
#'
#' @return TRUE or FALSE
#' @concept helper
different <- function(a, b) {
    !isTRUE(all.equal(a, b,
        check.attributes = FALSE, scale = 1,
        tolerance = 10 * .Machine$double.eps
    ))
}
