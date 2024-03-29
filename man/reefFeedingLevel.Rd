% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reef-project_methods.R
\name{reefFeedingLevel}
\alias{reefFeedingLevel}
\title{Reef feeding level}
\usage{
reefFeedingLevel(params, n, n_pp, n_other, t, encounter, ...)
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

\item{encounter}{A two dimensional array (predator species x predator size)
with the encounter rate.}

\item{...}{Unused}
}
\value{
A two dimensional array (predator species x predator size) with the
feeding level.
}
\description{
This function replaces the usual \code{\link[=mizerFeedingLevel]{mizerFeedingLevel()}} function and returns
the a feeding level of 0 for piscivores.
}
\seealso{
Other mizer rate functions: 
\code{\link{reefDegrade}()},
\code{\link{reefEncounter}()},
\code{\link{reefMort}()},
\code{\link{reefPredMort}()},
\code{\link{reefRates}()},
\code{\link{reefVulnerable}()}
}
\concept{extmort}
\concept{mizer rate functions}
