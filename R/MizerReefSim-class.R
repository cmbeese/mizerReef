#' Extension of MizerSim for reef-specific simulation results
#'
#' Stores simulation results and time-varying refuge profiles for reef models.
#'
#' @slot vulnerable List or array of refuge profiles for each time step
#' @seealso [MizerSim]
#' @export
#' @docType class
#' @name MizerReefSim
#' @concept simulation
#' @family simulation
setClass(
  "MizerReefSim",
  contains = "MizerSim",
  slots = c(
    vulnerable = "list" # can be matrix, vector, etc.
  )
)
