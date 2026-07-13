#' S4 marker class for mizerReef extension
#'
#' A formal S4 class that extends [MizerParams] but adds no new slots.
#' All reef-specific state (refuge, algae, detritus parameters) lives in
#' `other_params(params)`. The class label enables S3 dispatch of mizer
#' generics to mizerReef-specific methods.
#'
#' The class is **not** defined statically. Instead mizer creates it when the
#' package is loaded: `.onLoad()` calls [mizer::registerExtension()], which
#' recognises mizerReef as a dispatching extension from the S3 methods it
#' registers for its marker class and inserts `mizerReef` at the correct place
#' in the S4 hierarchy relative to any other extension packages loaded in the
#' same session. This lets mizerReef be chained with other extensions (for
#' example mizerMR, to add multiple resources) in either load order. A static
#' `contains = "MizerParams"` definition would fix mizerReef as a direct sibling
#' of every other extension and prevent such chaining, because a sealed class
#' cannot be re-parented.
#'
#' @seealso [MizerParams], [newReefParams()]
#' @docType class
#' @name mizerReef-class
#' @concept parameters
#' @family parameters
#' @keywords internal
NULL

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

    # Only unstructured components (algae, detritus) hold a single biomass
    # value per time step. Size-structured components contributed by other
    # extensions in the chain - such as the "MR" resource array added by
    # mizerMR - store a whole spectrum per time step and are reported by those
    # extensions instead, so they must be excluded here.
    scalar <- vapply(seq_len(ncol(sim@n_other)),
                     function(j) all(lengths(sim@n_other[, j]) == 1L),
                     logical(1))

    if (any(scalar)) {
        comp_mat <- matrix(unlist(sim@n_other[, scalar, drop = FALSE]),
            nrow = nrow(sim@n_other), ncol = sum(scalar)
        )
        dimnames(comp_mat) <- list(rownames(sim@n_other),
                                   colnames(sim@n_other)[scalar])

        b <- cbind(b, comp_mat[rownames(b), , drop = FALSE])
        names(dimnames(b)) <- dimname_names
    }

    mizer::ArrayTimeBySpecies(b,
        value_name = "Biomass", units = "g",
        params = sim@params
    )
}
