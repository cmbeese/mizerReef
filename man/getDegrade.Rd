% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reef-rates.R
\name{getDegrade}
\alias{getDegrade}
\title{Get vulnerability level at in time range t}
\usage{
getDegrade(object, n, n_pp, n_other, time_range, drop = TRUE, ...)
}
\arguments{
\item{object}{A \code{MizerParams} object or a \code{MizerSim} object}

\item{n}{A matrix of species abundances (species x size).}

\item{n_pp}{A vector of the resource abundance by size}

\item{n_other}{A list of abundances for other dynamical components of the
ecosystem}

\item{time_range}{A vector of times. Only the range of times is relevant,
i.e., all times between the smallest and largest will be selected.
The time_range can be character or numeric.}

\item{drop}{If \code{TRUE} then any dimension of length 1 will be removed
from the returned array.}

\item{...}{Unused}
}
\value{
If a \code{MizerParams} object is passed in, the function returns a
2 dimensional array (refuge size bin x refuge density) with the
new refuge densities for each size bin

\if{html}{\out{<div class="sourceCode">}}\preformatted{     If a `MizerSim` object is passed in, the function returns a three
     dimensional array (time step x refuge size bin x refuge density) 
     with the refuge density calculated at every time step in the 
     simulation. If \code{drop = TRUE} then the dimension of length 1 
     will be removed from the returned array.
}\if{html}{\out{</div>}}
}
\description{
Returns the proportion of fish at size \eqn{w} that are not hidden in
predation refuge and thus vulnerable to being encountered by predators.
}
\details{
This function uses \code{\link[=reefVulnerable]{reefVulnerable()}} to calculate the vulnerability to
predation.
}
\section{Setting the refuge profile}{


Refuge profiles account for the protective behavior of prey living in
high complexity environments (e.g. coral reefs) with access to predation
refuge. The refuge profile defines the proportion of fish within
user-defined length bins that are protected from being encountered
by a predator.

mizerReef determines how much fish weight a refuge can hold by converting
user provided fish length bins to weight bins with average length to weight
conversion parameters \eqn{\bar{a}} and \eqn{\bar{b}}, which describe
the length to weight relationship for the fish dummies used during data
collection. The default values are \eqn{\bar{a} = 0.025, \bar{b} = 3}.
Assuming that fish in shelters are neutrally buoyant, the mass of the
largest fish in each size category is proportional to the volume of
refuges classified
in that category.

A unique refuge profile is generated for each predator group x
prey group x prey size combination based on the given refuge profile
parameters as well as four values from \code{params@species_params}: length
to weight conversion values \code{a} and \code{b}, \code{refuge_user}, which is true
for groups utilize that predation refuge, and \code{bad_pred}, which is false
for predator groups whose body shape or predatory strategy allow them to
access fish within refuge (e.g. eels).

To ensure some food is always available to predators, the maximum
proportion of fish protected by refuge in any size class is set by
\code{max_protect}.

The refuge profile is used when calculating the food encounter rate in
\code{\link[=reefEncounter]{reefEncounter()}} and the predation mortality rate in \code{\link[=reefPredMort]{reefPredMort()}}.
Its entries are dimensionless values between 0 and 1 which represent the
proportion of fish in the corresponding prey and size categories that are
hidden within refuge and thus cannot be encountered by predators. If no
refuge is available then predator-prey interactions are determined
entirely by size-preference.

The mizerReef package provides three methods to define the refuge profile.

\itemize{

\item \strong{Sigmoidal Method}: \cr
This method is preferred for data-poor reefs or reefs
where the refuge distribution is unknown. It is also ideal for systems
where only one functional group is expected to be utilizing refuge. The
proportion of fish with access to refuge \eqn{ R_j(w_p) } is given by

\deqn{ R_j(w_p) = \frac{r}
                        {1 + e^{ \left( \Delta (w - W_{refuge}) \right)} }}{
      R_j(w_p) = r / (1 + e^{( Δ (w - W_{refuge}))})}

Here \eqn{W_{refuge}} marks the body weight at which refuge becomes
scarcer for prey. \eqn{r} defines the maximum proportion of fish with
access to predation refuge and is always less than or equal to
\code{max_protect}. \eqn{\alpha} controls the rate at which the
availability of refuge decreases with increasing body size. It
defaults to a steep slope of 100.

For this method, \code{method_params} should contain columns named
\code{prop_protect} and \code{L_refuge} that give the values for \eqn{r}
and the length at which refuge becomes scarce in cm.

\item \strong{Binned Method}: \cr
This method is appropriate for theoretical applications
and does not rely on empirical data. It sets refuge to a constant
proportion of fish within a given size range. The  proportion of fish
in group \eqn{j} with access to refuge is given by

\deqn{ R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~] }{
      R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~] }
   
where \eqn{r_k} is the proportion of fish with access to refuge in
size class \eqn{k}.

For this method, \code{method_params} should contain columns named
\code{start_L}and \code{end_L} which contain the starting and ending lengths \link{cm}
of each size bin and \code{prop_protect}, the proportion of fish protected
within each corresponding size bin.

\item \strong{Competitive Method}: \cr
This method is appropriate when refuge density
data is available for the modelled reef. The refuge density describes
the distribution of refuges \eqn{(no./m^2)} across defined fish body size
categories. The proportion of fish in size class \eqn{k} with access
to refuge is given by

\deqn{R_{j}(w_p) = \tau \cdot \frac{ \eta_{k} }
                                    { \sum_i \int_{w_{k-1}}^{w-k} N_i(w) \, dw}}{
      R_{j}(w_p) = \tau \eta_{k} / 
                  ( \sum_i \int_{w_{k-1}}^{w-k} N_i(w) \, dw ) }

where \eqn{ \tau } is the proportion of fish with access to refuge that
are expected to actually utilize  it, \eqn{ \eta_{k}} is the density of
refuges in size range \eqn{(w_{k-1}, w_k]} and
\eqn{\sum_{i} \int_{w_{k-1}}^{w_k} N_i(w)~dw} gives the density
of fish from any group in size range \eqn{(w_{k-1}, w_k]}.
This represents the density of competitors for refuges in
size class \eqn{k}.

For this method, \code{method_params} should contain columns named
\code{start_L}and \code{end_L} which contain the starting and ending lengths \link{cm}
of each size bin and \code{refuge_density}, the number of refuges available
in each size bin (no/m^2).

}

Users can also set a noncomplex reef with no habitat refuge. This option is
convenient for finding steady state parameters.

This function checks that the supplied refuge parameters are valid, adds
relevant columns to the \code{species_params} data frame, and stores refuge
parameters in the \code{other_params} slot of the \code{params} object.

Refuge profile parameters can be input in a spreadsheet program and saved
as a .csv file. The data can then be read into R using the command
\code{read.csv()}.
}

\seealso{
Other rate functions: 
\code{\link{getSenMort}()},
\code{\link{getVulnerable}()}
}
\concept{degradation}
\concept{rate functions}
