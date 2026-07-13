#' mizerReef: A mizer extension package for modelling structurally complex marine ecosystems.
#' Captures the role of habitat complexity as predation refuge and add benthic coupling.
#'
#' @description This is an extension package for the mizer package
#' (https://sizespectrum.org/mizer/) that makes it easy to set up a mizer model
#' with predation refuges, detritus, and algae, thereby allowing for more
#' realistic modelling of coral reef fisheries.
#'
#'This package was developed to support the creation and exploration of a model
#'for coral reefs where we could compare areas based on their benthic complexity.
#'
#' @import mizer mizerExperimental ggplot2
#' @importFrom plyr aaply
#' @importFrom plotly ggplotly
#' @importFrom lubridate now
#' @importFrom methods is
#' @importFrom assertthat assert_that is.flag is.number
#' @importFrom lifecycle deprecated
#' @importFrom dplyr %>% mutate left_join
#' @importFrom stats mvfft complete.cases
#' @md
#' @keywords internal
"_PACKAGE"

.onLoad <- function(libname, pkgname) {
    mizer::registerExtension(pkgname,
                             requirement = "sizespectrum/mizerReef")
    if (exists("caribbean_3_model", envir = asNamespace(pkgname), inherits = FALSE)) {
        ns <- asNamespace(pkgname)
        raw_caribbean_3_model <- get("caribbean_3_model", envir = ns)
        makeActiveBinding("caribbean_3_model",
                          fun = function() mizer::coerceToExtensionClass(raw_caribbean_3_model),
                          env = ns)
    }
}
