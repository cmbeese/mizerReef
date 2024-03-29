% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algae_components.R
\name{algae_dynamics}
\alias{algae_dynamics}
\title{Algae dynamics}
\usage{
algae_dynamics(params, n, n_other, rates, dt, ...)
}
\arguments{
\item{params}{A \link{MizerParams} object}

\item{n}{A matrix of current species abundances (species x size)}

\item{n_other}{Other dynamic components.}

\item{rates}{A list of rates as returned by \code{\link[=getRates]{getRates()}}}

\item{dt}{Time step size}

\item{...}{Unused}
}
\value{
A single number giving the algae biomass at next time step
}
\description{
Calculates the algal biomass at the next time step from the current
algae biomass
}
\details{
The time evolution of the algal biomass \eqn{B} is described by

\deqn{dB_A/dt = P_A - c_A \cdot B_A}{
      dB_A/dt = P_A - c_A * B_A}

where  \eqn{c_A} is the mass-specific rate of consumption calculated
with \code{algae_consumption()} and \eqn{P_A} is the rate at which algae
grows, calculated with \code{getAlgaeProduction()}.

The dynamical equation is solved analytically to

\deqn{B_A(t+dt) = B_A(t) e^{(- c_A \cdot dt)}
             +\frac{P_A}{c_A}
             (1- e^{(-c_A \cdot dt)}).}{
      B_A(t+dt) = B(t) exp(-c_A * dt) 
              + p_A/c_A * (1 - exp(-c_A * dt)).}

This avoids the stability problems that would arise if we used the Euler
method to solve the equation numerically.
}
\seealso{
\code{\link[=detritus_dynamics]{detritus_dynamics()}}, \code{\link[=algae_consumption]{algae_consumption()}},
\code{\link[=getAlgaeConsumption]{getAlgaeConsumption()}}, \code{\link[=getAlgaeProduction]{getAlgaeProduction()}}
}
\concept{algae}
