% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reef-summary_plots.R
\name{plotProductivity}
\alias{plotProductivity}
\alias{plotlyProductivity}
\title{Plot the total productivity for each functional group}
\usage{
plotProductivity(
  object,
  species = NULL,
  return_data = FALSE,
  min_fishing_l = NULL,
  max_fishing_l = NULL,
  start_time = NULL,
  end_time = NULL,
  ...
)

plotlyProductivity(object, species = NULL, ...)
}
\arguments{
\item{object}{An object of class \linkS4class{MizerParams}}

\item{species}{The species to be selected. Optional. By default all
species are selected. A vector of species names, or a
numeric vector with the species indices, or a logical
vector indicating for each species whether it is to be
selected (TRUE) or not.}

\item{return_data}{A boolean value that determines whether the formatted
data used for the plot is returned instead of the plot
itself. Default value is FALSE.}

\item{min_fishing_l}{parameters be passed to \code{\link[=getProductivity]{getProductivity()}}. The
minimum length (cm) of fished individuals for
productivity estimates. Defaults to 7 cm.}

\item{max_fishing_l}{parameters be passed to \code{\link[=getProductivity]{getProductivity()}}. The
maximum length (cm) of fished individuals for
productivity estimates. Defaults to max length.}

\item{start_time}{The first time to be plotted. Default is the beginning
of the time series.}

\item{end_time}{The last time to be plotted. Default is the end of the time
series.}

\item{...}{unused}
}
\value{
A ggplot2 object, unless \code{return_data = TRUE}, in which case a data
frame with the the productivity for each functional group
is returned.
}
\description{
When called with a \linkS4class{MizerParams} object the steady state
productivity is plotted. When called with a \linkS4class{MizerSim}
object the productivity is plotted.
}
\seealso{
\code{\link[=plotBiomass]{plotBiomass()}}, \code{\link[=plot2TotalBiomass]{plot2TotalBiomass()}}, \code{\link[=plotTotalBiomassRelative]{plotTotalBiomassRelative()}},
\code{\link[=plotProductivity]{plotProductivity()}}, \code{\link[=plot2Productivity]{plot2Productivity()}}, \code{\link[=plotProductivityRelative]{plotProductivityRelative()}}

Other plotting functions: 
\code{\link{plot2Productivity}()},
\code{\link{plot2TotalBiomass}()},
\code{\link{plotBiomass}()},
\code{\link{plotProductivityRelative}()},
\code{\link{plotRefuge}()},
\code{\link{plotTotalBiomassRelative}()},
\code{\link{plotTotalBiomass}()},
\code{\link{plotVulnerable}()}
}
\concept{plotting functions}
\concept{sumplots}
