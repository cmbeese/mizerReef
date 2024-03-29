% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reef-project_methods.R
\name{reefDegrade}
\alias{reefDegrade}
\title{Scales the refuge density by a given value at set times}
\usage{
reefDegrade(params, n, n_pp, n_other, t, ...)
}
\arguments{
\item{params}{A \linkS4class{MizerParams} object}

\item{n}{A matrix of species abundances (species x size).}

\item{n_pp}{A vector of the resource abundance by size}

\item{n_other}{A list of abundances for other dynamical components of the
ecosystem}

\item{t}{The time for which to do the calculation (Not used by standard
mizer rate functions but useful for extensions with time-dependent
parameters.)}

\item{...}{Unused}
}
\value{
A new methods parameters data frame scaled by bleaching
}
\description{
Allows for the degradation of coral reef habitat structure following an
acute disturbance by decreasing the availability of refuge over time.
}
\seealso{
Other mizer rate functions: 
\code{\link{reefEncounter}()},
\code{\link{reefFeedingLevel}()},
\code{\link{reefMort}()},
\code{\link{reefPredMort}()},
\code{\link{reefRates}()},
\code{\link{reefVulnerable}()}
}
\concept{degradation}
\concept{mizer rate functions}
