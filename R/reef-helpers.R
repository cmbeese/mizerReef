#' Remove some species from the model
#'
#' This calls `mizer::removeSpecies()` and in addition removes the relevant
#' row from the detritus and algae consumption arrays `rho`.
#' @param params A MizerParams object
#' @param species The species to be removed. A vector of species names, or a
#'   numeric vector of species indices, or a logical vector indicating for each
#'   species whether it is to be removed (TRUE) or not.
#' @return A MizerParams object with fewer species.
#' @examples
#' params <- NWMed_params
#' species_params(params)$species
#' params <- removeSpecies(params, c("Poor cod", "Horse mackerel"))
#' species_params(params)$species
#' @export
removeSpecies <- function(params, species) {
    p <- mizer::removeSpecies(params, species)
    species <- valid_species_arg(params, species,
                                 return.logical = TRUE)
    keep <- !species

    # Remove algae rho values for species
    p@other_params$algae$rho <-
        p@other_params$algae$rho[keep, , drop = FALSE]

    # Remove detritus rho values for species
    p@other_params$detritus$rho <-
        p@other_params$detritus$rho[keep, , drop = FALSE]

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
    !isTRUE(all.equal(a, b, check.attributes = FALSE, scale = 1,
                      tolerance = 10 * .Machine$double.eps))
}
