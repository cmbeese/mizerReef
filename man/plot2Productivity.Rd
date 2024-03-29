% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reef-summary_plots.R
\name{plot2Productivity}
\alias{plot2Productivity}
\alias{plotly2Productivity}
\title{Plot the fisheries productivity of two models or two different size ranges
in the same plot}
\usage{
plot2Productivity(
  object1,
  object2,
  species = NULL,
  name1 = "First",
  name2 = "Second",
  min_fishing_l1 = NULL,
  max_fishing_l1 = NULL,
  min_fishing_l2 = NULL,
  max_fishing_l2 = NULL,
  stack = FALSE,
  return_data = FALSE,
  ...
)

plotly2Productivity(object1, object2, ...)
}
\arguments{
\item{object1}{First MizerParams or MizerSim object.}

\item{object2}{Second MizerParams or MizerSim object.}

\item{species}{The groups to be selected. Optional. By default all target
groups are selected. A vector of groups names, or a numeric
vector with the groups indices, or a logical vector
indicating for each group whether it is to be selected
(TRUE) or not.}

\item{name1}{An optional string with the name for the first model, to be
used in the legend. Set to "First" by default.}

\item{name2}{An optional string with the name for the second model, to be
used in the legend. Set to "Second" by default.}

\item{min_fishing_l1}{Optional.  The minimum length (cm) of fished
individuals for model 2. Defaults to 7cm.
A parameter passed to \code{\link[=getProductivity]{getProductivity()}}.}

\item{max_fishing_l1}{Optional.  The maximum length (cm) of fished
individuals for model 1. Defaults to max length.
A parameter passed to \code{\link[=getProductivity]{getProductivity()}}.}

\item{min_fishing_l2}{Optional.  The minimum length (cm) of fished
individuals for model 2. Defaults to 7cm.
A parameter passed to \code{\link[=getProductivity]{getProductivity()}}.}

\item{max_fishing_l2}{Optional.  The maximum length (cm) of fished
individuals for model 1. Defaults to max length.
A parameter passed to \code{\link[=getProductivity]{getProductivity()}}.}

\item{stack}{A boolean value that determines whether bars are separated
by species. Defaults to FALSE. If true, returns a stacked
barplot with the total biomass for each group instead of
individual bars for each group. Useful for comparison
between steady states.}

\item{return_data}{A boolean value that determines whether the formatted
data used for the plot is returned instead of the plot
itself. Default value is FALSE.}

\item{...}{
  Arguments passed on to \code{\link[=plotProductivity]{plotProductivity}}
  \describe{
    \item{\code{object}}{An object of class \linkS4class{MizerParams} or
\linkS4class{MizerSim}}
    \item{\code{start_time}}{The first time to be plotted. Default is the beginning
of the time series.}
    \item{\code{end_time}}{The last time to be plotted. Default is the end of the time
series.}
    \item{\code{facet}}{A boolean value indicating whether to facet the result plot
by species group. Defaults to TRUE.}
    \item{\code{total}}{A boolean value that determines whether the total productivity
from all species is plotted as well. Default is TRUE.}
    \item{\code{min_fishing_l}}{parameters be passed to \code{\link[=getProductivity]{getProductivity()}}. The
minimum length (cm) of fished individuals for
productivity estimates. Defaults to 7 cm.}
    \item{\code{max_fishing_l}}{parameters be passed to \code{\link[=getProductivity]{getProductivity()}}. The
maximum length (cm) of fished individuals for
productivity estimates. Defaults to max length.}
  }}
}
\value{
A ggplot2 object, unless \code{return_data = TRUE}, in which case a data
frame with the the productivity for each Species Group
by model is returned.
}
\description{
When called with a \linkS4class{MizerParams}
object the steady state productivities are plotted.
}
\seealso{
\code{\link[=plotBiomass]{plotBiomass()}}, \code{\link[=plot2TotalBiomass]{plot2TotalBiomass()}}, \code{\link[=plotTotalBiomassRelative]{plotTotalBiomassRelative()}},
\code{\link[=plotProductivity]{plotProductivity()}}, \code{\link[=plot2Productivity]{plot2Productivity()}}, \code{\link[=plotProductivityRelative]{plotProductivityRelative()}}

Other plotting functions: 
\code{\link{plot2TotalBiomass}()},
\code{\link{plotBiomass}()},
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
