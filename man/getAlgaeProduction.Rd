% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algae_components.R
\name{getAlgaeProduction}
\alias{getAlgaeProduction}
\title{Algae production rate}
\usage{
getAlgaeProduction(params)
}
\arguments{
\item{params}{MizerParams}
}
\value{
The annual growth rate of algae per square meter
}
\description{
This is the rate in grams/year/m^-2 at which the system produces algae
biomass. The rate is set so that production and consumption are equal for
chosen steady state abundances.
}
\seealso{
\code{\link[=getAlgaeConsumption]{getAlgaeConsumption()}}, \code{\link[=algae_dynamics]{algae_dynamics()}}, \code{\link[=getAlgaeProduction]{getAlgaeProduction()}}
}
\concept{algae}
