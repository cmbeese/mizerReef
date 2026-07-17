#' Upgrade a mizerReef params object to the current layout
#'
#' Called automatically by [mizer::validParams()] when an object created with
#' an older version of mizerReef is loaded. Migrates data stored in the old
#' custom S4 slots (`refuge_params`, `algae_params`, `detritus_params`) to the
#' `other_params` sub-lists used since version 2.1.0.
#'
#' @param object A `mizerReef` params object (possibly from an older version).
#' @param ... Unused.
#' @return An upgraded `mizerReef` object.
#' @importFrom methods .hasSlot slot
#' @exportS3Method utils::upgrade
upgrade.mizerReef <- function(object, ...) {
    # Migrate old S4 slot data to other_params sub-lists if needed.
    # The old class had three extra slots: refuge_params, algae_params,
    # and detritus_params. In mizerReef >= 2.1.0 these are stored as named
    # lists inside other_params instead.
    #
    # Because the slots no longer exist in the class definition, attempting to
    # access them via @ would raise an error. We use tryCatch to detect the
    # old layout structurally and migrate it safely.
    tryCatch({
        old_refuge <- .hasSlot(object, "refuge_params")
        old_algae  <- .hasSlot(object, "algae_params")
        old_det    <- .hasSlot(object, "detritus_params")
        if (old_refuge && is.null(object@other_params$refuge_params)) {
            object@other_params$refuge_params  <- slot(object, "refuge_params")
        }
        if (old_algae && is.null(object@other_params$algae_params)) {
            object@other_params$algae_params   <- slot(object, "algae_params")
        }
        if (old_det && is.null(object@other_params$detritus_params)) {
            object@other_params$detritus_params <- slot(object, "detritus_params")
        }
    }, error = function(e) NULL)  # silently ignore if slots do not exist
    object
}
