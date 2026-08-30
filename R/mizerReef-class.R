#' mizerReef extension classes
#'
#' S3 extension classes for [MizerParams] and [MizerSim] that enable S3 dispatch
#' for extension-specific methods.
#'
#' The class names are ordinary entries in the object's S3 class vector. All
#' reef-specific data lives in `other_params(params)` or in component
#' parameters (see [setComponent()]).
#'
#' Objects of class `mizerReef` are created by [newReefParams()].
#' Objects of class `mizerReefSim` are returned automatically by [project()]
#' when called on a `mizerReef` params object.
#'
#' No class declaration is needed. [newReefParams()] records the
#' extension on the object with [mizer::recordExtension()] and then calls
#' [mizer::coerceToExtensionClass()].
#'
#' @seealso [MizerParams], [newReefParams()]
#' @name mizerReef-class
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
