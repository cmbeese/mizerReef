% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reef-summary_plots.R
\name{plotBiomass}
\alias{plotBiomass}
\alias{plotlyBiomass}
\title{Plot the biomass of Species Groups and unstructured components through
time}
\usage{
plotBiomass(
  sim,
  species = NULL,
  start_time,
  end_time,
  y_ticks = 6,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  highlight = NULL,
  return_data = FALSE,
  min_length = NULL,
  ...
)

plotlyBiomass(
  sim,
  species = NULL,
  start_time,
  end_time,
  y_ticks = 6,
  ylim = c(NA, NA),
  total = FALSE,
  background = TRUE,
  highlight = NULL,
  min_length = NULL,
  ...
)
}
\arguments{
\item{sim}{An object of class \linkS4class{MizerSim}}

\item{species}{The groups to be selected. Optional. By default all target
groups are selected. A vector of groups names, or a numeric
vector with the groups indices, or a logical vector
indicating for each group whether it is to be selected
(TRUE) or not.}

\item{start_time}{The first time to be plotted. Default is the beginning
of the time series.}

\item{end_time}{The last time to be plotted. Default is the end of the time
series.}

\item{y_ticks}{The approximate number of ticks desired on the y axis.}

\item{ylim}{A numeric vector of length two providing lower and upper limits
for the y axis. Use NA to refer to the existing minimum or
maximum. Any values below 1e-20 are always cut off.}

\item{total}{A boolean value that determines whether the total biomass from
species is plotted as well. Default is FALSE.}

\item{background}{A boolean value that determines whether background
species are included. Ignored if the model does not
contain background species. Default is TRUE.}

\item{highlight}{Name or vector of names of the species to be highlighted.}

\item{return_data}{A boolean value that determines whether the formatted
data used for the plot is returned instead of the plot
itself. Default value is FALSE}

\item{min_length}{minimum length of organism to include in biomass
calculation}

\item{...}{
  Arguments passed on to \code{\link[mizer:get_size_range_array]{mizer::get_size_range_array}}
  \describe{
    \item{\code{min_w}}{Smallest weight in size range. Defaults to smallest weight in
the model.}
    \item{\code{max_w}}{Largest weight in size range. Defaults to largest weight in the
model.}
    \item{\code{min_l}}{Smallest length in size range. If supplied, this takes
precedence over \code{min_w}.}
    \item{\code{max_l}}{Largest length in size range. If supplied, this takes precedence
over \code{max_w}.}
  }}
}
\value{
A ggplot2 object, unless \code{return_data = TRUE}, in which case a data
frame with the four variables 'Year', 'Biomass', 'Species',
'Legend' is returned.
}
\description{
After running a projection, the biomass of each Species Group and each
unstructured component can be plotted against time. The biomass is
calculated within user defined size limits (min_w, max_w, min_l, max_l, see
\code{\link[=get_size_range_array]{get_size_range_array()}}).
}
\seealso{
\code{\link[=plotBiomass]{plotBiomass()}}, \code{\link[=plot2TotalBiomass]{plot2TotalBiomass()}}, \code{\link[=plotTotalBiomassRelative]{plotTotalBiomassRelative()}},
\code{\link[=plotProductivity]{plotProductivity()}}, \code{\link[=plot2Productivity]{plot2Productivity()}}, \code{\link[=plotProductivityRelative]{plotProductivityRelative()}}

Other plotting functions: 
\code{\link{plot2Productivity}()},
\code{\link{plot2TotalBiomass}()},
\code{\link{plotProductivity}()},
\code{\link{plotProductivityRelative}()},
\code{\link{plotRefuge}()},
\code{\link{plotRelativeContribution}()},
\code{\link{plotSpectraRelative}()},
\code{\link{plotTotalAbundance}()},
\code{\link{plotTotalBiomass}()},
\code{\link{plotTotalBiomassRelative}()},
\code{\link{plotVulnerable}()}
}
\concept{plotting functions}
\concept{sumplots}
