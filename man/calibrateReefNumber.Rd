% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reef-components.R
\name{calibrateReefNumber}
\alias{calibrateReefNumber}
\title{Calibrate the model scale to match total observed number}
\usage{
calibrateReefNumber(params)
}
\arguments{
\item{params}{A MizerParams object}
}
\value{
A MizerParams object
}
\description{
Replaces mizer's \code{\link[mizer:calibrateNumber]{mizer::calibrateNumber()}} function. Given a MizerParams
object \code{params} for which number observations are available for at least
some species via the \code{number_observed} column in the species_params data
frame, this function returns an updated MizerParams object which is
rescaled with \code{\link[=scaleReefModel]{scaleReefModel()}} so that the total number in
the model agrees with the total observed number.
}
\details{
Number observations usually only include individuals above a certain size.
This size should be specified in a number_cutoff column of the species
parameter data frame. If this is missing, it is assumed that all sizes are
included in the observed number, i.e., it includes larval number.

After using this function the total number in the model will match the
total number, summed over all species. However the numbers of the
individual species will not match observations yet, with some species
having numbers that are too high and others too low. So after this
function you may want to use \code{\link[=matchNumbers]{matchNumbers()}}. This is described in the
blog post at https://bit.ly/2YqXESV.

If you have observations of the yearly yield instead of numbers, you can
use \code{\link[=calibrateYield]{calibrateYield()}} instead of this function.
}
\concept{calibration}
