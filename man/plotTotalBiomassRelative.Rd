% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reef-summary_plots.R
\name{plotTotalBiomassRelative}
\alias{plotTotalBiomassRelative}
\alias{plotlyTotalBiomassRelative}
\title{Plot the relative difference in between the total fishable biomasses of each
each Species Group at steady state}
\usage{
plotTotalBiomassRelative(
  object1,
  object2,
  diff_method,
  species = NULL,
  min_fishing_l1 = NULL,
  max_fishing_l1 = NULL,
  min_fishing_l2 = NULL,
  max_fishing_l2 = NULL,
  return_data = FALSE,
  ...
)

plotlyTotalBiomassRelative(object1, object2, diff_method, ...)
}
\arguments{
\item{object1}{First MizerParams or MizerSim object.}

\item{object2}{Second MizerParams or MizerSim object.}

\item{diff_method}{The method to calculate the relative change between
models. If \code{percent.change}, the percent change is
calculated relative to the value from object 1 with
formula 100*(new-old)/old. If \code{rel.diff} the relative
difference is returned given by (new - old)/(old + new).}

\item{species}{The groups to be selected. Optional. By default all target
groups are selected. A vector of groups names, or a numeric
vector with the groups indices, or a logical vector
indicating for each group whether it is to be selected
(TRUE) or not.}

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

\item{return_data}{A boolean value that determines whether the formatted
data used for the plot is returned instead of the plot
itself. Default value is FALSE.}

\item{...}{
  Arguments passed on to \code{\link[=plotTotalBiomass]{plotTotalBiomass}}
  \describe{
    \item{\code{object}}{An object of class \linkS4class{MizerParams}}
    \item{\code{min_fishing_l}}{The minimum length (cm) for biomass estimates.
Defaults to smallest size.}
    \item{\code{max_fishing_l}}{The maximum length (cm) of for biomass estimates.
Defaults to max length.}
  }}
}
\value{
A ggplot2 object, unless \code{return_data = TRUE}, in which case a data
frame with the the total steady state biomass for each functional
group by model is returned as well as another column called
\code{rel_diff}that gives the relative difference between the two
values.
}
\description{
This functions creates a barplot with the relative change in biomass of
each Species Group within a size range between either (1) two different
mizerParams objects (two models) or (2) two different size ranges.
}
\details{
This function is usually used in conjunction with
\code{\link[=plotProductivityRelative]{plotProductivityRelative()}} to check for decoupling between biomass and
productivity.

The individual productivity rates are calculated by the
\code{\link[=plotTotalBiomass]{plotTotalBiomass()}} function which is passed all additional arguments you
supply. See \code{\link[=plotTotalBiomass]{plotTotalBiomass()}} for more details.

To compare between different size ranges, use the \code{min_fishing_l1}
and \code{max_fishing_l1} arguments for the first size range and  the
\code{min_fishing_l2}and \code{max_fishing_l2} arguments for the second.
}
\seealso{
\code{\link[=plotBiomass]{plotBiomass()}}, \code{\link[=plot2TotalBiomass]{plot2TotalBiomass()}}, \code{\link[=plotTotalBiomassRelative]{plotTotalBiomassRelative()}},
\code{\link[=plotProductivity]{plotProductivity()}}, \code{\link[=plot2Productivity]{plot2Productivity()}}, \code{\link[=plotProductivityRelative]{plotProductivityRelative()}}

Other plotting functions: 
\code{\link{plot2Productivity}()},
\code{\link{plot2TotalBiomass}()},
\code{\link{plotBiomass}()},
\code{\link{plotProductivity}()},
\code{\link{plotProductivityRelative}()},
\code{\link{plotRefuge}()},
\code{\link{plotRelativeContribution}()},
\code{\link{plotSpectraRelative}()},
\code{\link{plotTotalAbundance}()},
\code{\link{plotTotalBiomass}()},
\code{\link{plotVulnerable}()}
}
\concept{plotting functions}
\concept{sumplots}
