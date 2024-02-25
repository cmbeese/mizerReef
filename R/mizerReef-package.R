#' mizerReef: A mizer extension package for modelling tropical coral reefs.
#' Includes benthic complexity through the provision of predation refuge.
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

