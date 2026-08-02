#' @title Deprecated: Use `newReefParams()` instead
#'
#' @description
#' `MizerReefParams()` is a deprecated constructor. Use [newReefParams()] to
#' create a mizerReef params object.
#'
#' @param ... Arguments passed to [MizerParams()]
#' @param refuge_params Ignored (now stored in `other_params`)
#' @param algae_params Ignored (now stored in `other_params`)
#' @param detritus_params Ignored (now stored in `other_params`)
#' @return A `mizerReef` object
#' @seealso [newReefParams()]
#' @export
#' @importFrom lifecycle deprecated
MizerReefParams <- function(...,
                            refuge_params = list(),
                            algae_params = list(),
                            detritus_params = list()) {
    lifecycle::deprecate_warn("2.0.0", "MizerReefParams()", "newReefParams()")
    newReefParams(...)
}
