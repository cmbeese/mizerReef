% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reef-project_methods.R
\name{reefEncounter}
\alias{reefEncounter}
\title{Get encounter rate needed to project a mizerReef model}
\usage{
reefEncounter(
  params,
  n,
  n_pp,
  n_other,
  t,
  vulnerable = reefVulnerable(params, n, n_pp, n_other, t, new_rd = reefDegrade(params,
    n, n_pp, n_other, t)),
  ...
)
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

\item{vulnerable}{A two dimensional array (prey species x prey size) with
the proportion of prey vulnerable to being encountered.}

\item{...}{Unused}
}
\value{
A named two dimensional array (predator species x predator size) with
the encounter rates.
}
\description{
Calculates the rate \eqn{E_i(w)} at which a predator from group \eqn{i} and
weight \eqn{w} encounters food (grams/year). You would not usually call this
function directly but instead use \code{\link[=getEncounter]{getEncounter()}}, which then calls this
function.
}
\section{Predation encounter}{


The encounter rate \eqn{E_i(w)} at which a predator of species \eqn{i}
and weight \eqn{w} encounters food has contributions from the encounter
of fish prey and of resources. This is determined by summing over all
prey species and the resource spectrum and then integrating over all
prey sizes \eqn{w_p}, weighted by predation kernel \eqn{\phi(w,w_p)}:

\deqn{
 E_i(w) = \gamma_i(w) \int
 \left( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij}
 V_{ji}(w_p) N_j(w_p) \right)
 \phi_i(w,w_p) w_p \, dw_p.
 }{\gamma_i(w) \int
 ( \theta_{ip} N_R(w_p) + \sum_{j} \theta_{ij} V_{ji}(w_p) N_j(w_p) )
 \phi_i(w,w_p) w_p dw_p.}
Here \eqn{N_j(w)} is the abundance density of species \eqn{j} and
\eqn{N_R(w)} is the abundance density of resource. The overall
prefactor \eqn{\gamma_i(w)} determines the predation power of the
predator. It could be interpreted as a search volume and is set with
the \code{\link[=setSearchVolume]{setSearchVolume()}} function.

The predation kernel \eqn{\phi(w,w_p)}is set with the \code{\link[=setPredKernel]{setPredKernel()}}
function.

The vulnerability to predation, \eqn{V_{ji}(w)} accounts for protective
behavior of the prey. The parameters that control this are set with the
\code{\link[=setRefuge]{setRefuge()}} function.

The species interaction matrix \eqn{\theta_{ij}} is set with
\code{\link[=setInteraction]{setInteraction()}} and the resource interaction vector \eqn{\theta_{ip}}
is taken from the \code{interaction_resource}column in
\code{params@species_params}.
}

\section{Details}{

The encounter rate is multiplied by \eqn{1-f_0} to obtain the consumption
rate, where \eqn{f_0} is the feeding level calculated with
\code{\link[=getFeedingLevel]{getFeedingLevel()}}. This is used by the \code{\link[=project]{project()}} function for performing
simulations.

The function returns values also for sizes outside the size-range of the
species. These values should not be used, as they are meaningless.

If your model contains additional components that you added with
\code{\link[=setComponent]{setComponent()}} and for which you specified an \code{encounter_fun} function then
the encounters of these components will be included in the returned value.
}

\seealso{
Other mizer rate functions: 
\code{\link{reefDegrade}()},
\code{\link{reefFeedingLevel}()},
\code{\link{reefMort}()},
\code{\link{reefPredMort}()},
\code{\link{reefRates}()},
\code{\link{reefVulnerable}()}
}
\concept{mizer rate functions}
\concept{refugeRates}
