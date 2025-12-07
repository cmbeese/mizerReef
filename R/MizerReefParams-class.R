#' Extension of MizerParams for reef-specific model parameters
#'
#' Stores reef-specific parameters, including predation refuge,
#' and algae and detritus resources.
#'
#' @slot refuge_params List for predation refuge parameters
#' @slot algae_params List for algae resource parameters
#' @slot detritus_params List for detritus resource parameters
#' @seealso [MizerParams]
#' @examples
#' # Example usage:
#' data(karpata_species)
#' params <- MizerReefParams(
#'   species_params = karpata_species,
#'   refuge_params = list(),
#'   algae_params = list(),
#'   detritus_params = list()
#' )
#' @export
#' @docType class
#' @name MizerReefParams
#' @concept parameters
#' @family parameters
#' @importFrom methods new
setClass(
  "MizerReefParams",
  contains = "MizerParams",
  slots = c(
    refuge_params = "list",
    algae_params = "list",
    detritus_params = "list"
  )
)

#' Constructor for MizerReefParams
#'
#' @description
#' Creates a `MizerReefParams` object with slots for reef-specific parameters,
#' including predation refuge and algae and detritus resources.
#'
#' @param ... Arguments for `MizerParams`
#' @param refuge_params List for refuge parameters
#' @param algae_params List for algae parameters
#' @param detritus_params List for detritus parameters
#' @return A `MizerReefParams` object
#' @seealso [MizerParams]
#' @export
MizerReefParams <- function(...,
                            refuge_params = list(),
                            algae_params = list(),
                            detritus_params = list()) {
  params <- MizerParams(...)
  new("MizerReefParams",
    params,
    refuge_params = refuge_params,
    algae_params = algae_params,
    detritus_params = detritus_params
  )
}
