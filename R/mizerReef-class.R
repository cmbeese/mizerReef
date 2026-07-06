#' S4 marker class for mizerReef extension
#'
#' A formal S4 class that extends [MizerParams] but adds no new slots.
#' All reef-specific state (refuge, algae, detritus parameters) lives in
#' `other_params(params)`. The class label enables S3 dispatch of mizer
#' generics to mizerReef-specific methods.
#'
#' @seealso [MizerParams], [newReefParams()]
#' @export
#' @docType class
#' @name mizerReef-class
#' @concept parameters
#' @family parameters
#' @importFrom methods setClass
setClass("mizerReef", contains = "MizerParams")

#' @rdname mizerReef-class
#' @export
setClass("mizerReefSim", contains = "MizerSim")

#' Get the biomass of species and unstructured components through time
#'
#' Extends [mizer::getBiomass()] to also include the algae and detritus
#' biomasses alongside the species biomasses, so that [mizer::plotBiomass()]
#' (which calls `getBiomass()` internally) shows them without needing its
#' own override.
#'
#' @param object A `mizerReefSim` object
#' @param ... Passed on to [mizer::getBiomass()]
#' @return An [mizer::ArrayTimeBySpecies] object (time x species/component)
#' @method getBiomass mizerReefSim
#' @export
getBiomass.mizerReefSim <- function(object, ...) {
    sim <- object
    b <- unclass(NextMethod())
    dimname_names <- names(dimnames(b))

    comp_mat <- matrix(unlist(sim@n_other),
        nrow = nrow(sim@n_other), ncol = ncol(sim@n_other)
    )
    dimnames(comp_mat) <- dimnames(sim@n_other)

    b <- cbind(b, comp_mat[rownames(b), , drop = FALSE])
    names(dimnames(b)) <- dimname_names

    mizer::ArrayTimeBySpecies(b,
        value_name = "Biomass", units = "g",
        params = sim@params
    )
}
