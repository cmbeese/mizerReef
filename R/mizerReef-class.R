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
