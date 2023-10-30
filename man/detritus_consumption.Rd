% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detritus_components.R
\name{detritus_consumption}
\alias{detritus_consumption}
\title{Mass-specific detritus consumption rate}
\usage{
detritus_consumption(params, n = params@initial_n, rates = getRates(params))
}
\arguments{
\item{params}{MizerParams}

\item{n}{A matrix of current species abundances (species x size)}

\item{rates}{A list of rates as returned by \code{\link[=getRates]{getRates()}}}
}
\value{
The mass-specific consumption rate of detritus in grams per year.
}
\description{
This mass-specific consumption rate is used in \code{detritus_dynamics()} to
calculate the detritus biomass at the next time step. To get the
non-mass-specific consumption rate, use \code{getDetritusConsumption()}.
}
\details{
The consumption rate by detritivorous fish is determined by
\code{other_params(params)$detritus$rho}
}