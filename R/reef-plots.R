library(ggplot2)
library(plotly)

#' Plot the biomass of functional groups and unstructured components through
#' time
#'
#' After running a projection, the biomass of each species and each
#' unstructured component can be plotted against time. The biomass is
#' calculated within user defined size limits (min_w, max_w, min_l, max_l, see
#' [get_size_range_array()]).
#'
#' @param sim An object of class \linkS4class{MizerSim}
#' @param species The groups to be selected. Optional. By default all target
#'   groups are selected. A vector of groups names, or a numeric vector with
#'   the groups indices, or a logical vector indicating for each groups
#'   whether it is to be selected (TRUE) or not.
#' @param y_ticks The approximate number of ticks desired on the y axis.
#' @param start_time The first time to be plotted. Default is the beginning of
#'   the time series.
#' @param end_time The last time to be plotted. Default is the end of the time
#'   series.
#' @param ylim A numeric vector of length two providing lower and upper limits
#'   for the y axis. Use NA to refer to the existing minimum or maximum. Any
#'   values below 1e-20 are always cut off.
#' @param total A boolean value that determines whether the total biomass from
#'   all species is plotted as well. Default is FALSE.
#' @param background A boolean value that determines whether background species
#'   are included. Ignored if the model does not contain background species.
#'   Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' @param return_data A boolean value that determines whether the formatted data
#'   used for the plot is returned instead of the plot itself. Default value is
#'   FALSE
#' @inheritDotParams mizer::get_size_range_array -params
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame with the four variables 'Year', 'Biomass', 'Species', 'Legend' is
#'   returned.
#'
#' @import ggplot2
#' @export
#' @concept plots
#' @family plotting functions
plotBiomass <- function(sim, species = NULL,
                        start_time, end_time,
                        y_ticks = 6, ylim = c(NA, NA),
                        total = FALSE, background = TRUE,
                        highlight = NULL, return_data = FALSE,
                        ...) {

    # If there are no unstructured components then call the mizer function
    if (is.null(getComponent(params, "algae"))) {
        if (is.null(getComponent(params, "detritus"))) {
            return(mizer::plotBiomass(sim, species = species,
                                      start_time = start_time, 
                                      end_time = end_time,
                                      y_ticks = y_ticks, ylim = ylim,
                                      total = total, background = background,
                                      highlight = highlight,
                                      return_data = return_data, ...))
        }
    }

    # User mizer function to create dataframe
    df <- mizer::plotBiomass(sim, species = species,
                             start_time = start_time, end_time = end_time,
                             y_ticks = y_ticks, ylim = ylim,
                             total = total, background = background,
                             highlight = highlight,
                             return_data = TRUE, ...)

    params <- sim@params
    species <- valid_species_arg(sim, species)
    if (missing(start_time)) start_time <-
        as.numeric(dimnames(sim@n)[[1]][1])
    if (missing(end_time)) end_time <-
        as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]])
    if (start_time >= end_time) {
        stop("start_time must be less than end_time")
    }

    # Algae
    alg_biomass <- rowSums(sweep(sim@n_other$algae, 2, 
                                 params@dw_full * params@w_full, "*"))
    times <- as.numeric(names(alg_biomass))
    sel_times <- (times >= start_time) & (times <= end_time)
    br <- data.frame(Year = times[sel_times],
                     Biomass = alg_biomass[sel_times],
                     Species = "Algae")

    # Detritus
    d_biomass <- rowSums(sweep(sim@n_other$detritus, 2, 
                               params@dw_full * params@w_full, "*"))
    times <- as.numeric(names(d_biomass))
    sel_times <- (times >= start_time) & (times <= end_time)
    br <- data.frame(Year = times[sel_times],
                     Biomass = d_biomass[sel_times],
                     Species = "Detritus")

    # Implement ylim and a minimal cutoff
    min_value <- 1e-20
    br <- br[br$Biomass >= min_value &
                 (is.na(ylim[1]) | br$Biomass >= ylim[1]) &
                 (is.na(ylim[2]) | br$Biomass <= ylim[2]), c(1, 2, 3)]
    br$Legend <- br$Species

    # other components
    bc <- unlist(sim@n_other)
    dim(bc) <- dim(sim@n_other)
    dimnames(bc) <- dimnames(sim@n_other)
    times <- as.numeric(dimnames(bc)[[1]])
    bc <- bc[(times >= start_time) & (times <= end_time), , drop = FALSE]
    bc <- melt(bc)

    # Implement ylim and a minimal cutoff and bring columns in desired order
    min_value <- 1e-20
    bc <- bc[bc$value >= min_value &
                 (is.na(ylim[1]) | bc$value >= ylim[1]) &
                 (is.na(ylim[2]) | bc$value <= ylim[2]), c(1, 3, 2)]
    names(bc) <- c("Year", "Biomass", "Species")
    bc$Legend <- bc$Species

    plot_dat <- rbind(df, bc, br)
    if (return_data) return(plot_dat)

    plotDataFrame(plot_dat, params, xlab = "Year", ylab = "Biomass [g]",
                  ytrans = "log10",
                  y_ticks = y_ticks, highlight = highlight,
                  legend_var = "Legend")
}

#' @rdname plotBiomass
#' @export
plotlyBiomass <- function(sim,
                          species = NULL,
                          start_time,
                          end_time,
                          y_ticks = 6,
                          ylim = c(NA, NA),
                          total = FALSE,
                          background = TRUE,
                          highlight = NULL,
                          ...) {
    argg <- c(as.list(environment()), list(...))
    ggplotly(do.call("plotBiomass", argg),
             tooltip = c("Species", "Year", "Biomass"))
}

#' Plot the vulnerability to predation of species by size
#'
#' When called with a \linkS4class{MizerParams}
#' object the initial vulnerability is plotted.
#'
#' @param param An object of class \linkS4class{MizerParams}
#' @param species The species to be selected. Optional. By default all
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept plots
#' @seealso [plotting_functions], [setRefuge()]
plotVulnerable <- function(object,
                           species = NULL,
                           ...) {

    if (is(object, "MizerSim")) {

        stop('This functionality is not set up yet you dumbass
             .')
    } else if (is(object, "MizerParams")) {

        params <- object
        vul <- getVulnerable(params)

        func_groups = c('Eels',
                        'Farming \n damsels',
                        'Herbivores',
                        'Parrotfish',
                        'Cryptobenthic \n predators',
                        'Engulfers',
                        'Grabbers',
                        'Nocturnal \n invertivores',
                        'Planktivores')

        # selector for desired species
        sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                    error_on_empty = TRUE)
        species <- dimnames(params@initial_n)$sp[sel_sp]
        species <- gsub('inverts',NA,species)
        species <- species[!is.na(species)]
        sel_sp <- which(!is.na(species))
        vul <- vul[sel_sp, , drop = FALSE]

        # Make data from from selected species
        plot_dat <- data.frame(w = rep(params@w, each = length(species)),
                               value = c(vul),
                               Species = species)

        # Set up colors
        legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
        plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
        linesize <- rep(0.8, length(legend_levels))
        names(linesize) <- names(params@linetype[legend_levels])

        p <- ggplot(plot_dat, aes(group = Species)) +
                    facet_wrap(~ Species, ncol = 3) +
                    theme(strip.background = element_blank(),
                          strip.text.x = element_blank())

        p + geom_line(aes(x = w, y = value,
                          colour = Legend, linetype = Legend,
                          linewidth = Legend)) +
            labs(colour = 'Functional Group', linetype = 'Functional Group',
                 linewidth = 'Functional Group') +
            scale_x_continuous(name = "Log Size [g]", trans = "log10",) +
            scale_y_continuous(name = "Proportion Protected", limits = c(0, 1)) +
            scale_colour_manual(values = params@linecolour[legend_levels],
                                labels = func_groups) +
            scale_linetype_manual(values = params@linetype[legend_levels],
                                  labels = func_groups) +
            scale_discrete_manual("linewidth", values = linesize,
                                  labels = c(func_groups))
    }
}


#' @rdname plotVulnerable
#' @export
plotlyVulnerable <- function(object,
                             species = NULL,
                             ...) {

    argg <- as.list(environment())
    ggplotly(do.call("plotVulnerable", argg),
             tooltip = c("Species", "w", "value"))
}

#' Plot the refuge profile, species by size
#'
#' When called with a \linkS4class{MizerParams}
#' object the initial vulnerability is plotted.
#'
#' @param param An object of class \linkS4class{MizerParams}
#' @param species The species to be selected. Optional. By default all
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept plots
#' @seealso [plotting_functions], [setRefuge()]
plotRefuge <- function(object, species = NULL,
                           ...) {

    if (is(object, "MizerSim")) {

        stop('This functionality is not set up yet you dumbass.')

    } else if (is(object, "MizerParams")) {

        params <- object

        vul <- getVulnerable(params)
        refuge <- (1-vul)

        # selector for desired species
        sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                    error_on_empty = TRUE)
        species <- dimnames(params@initial_n)$sp[sel_sp]
        species <- gsub('inverts',NA,species)
        species <- species[!is.na(species)]
        sel_sp <- which(!is.na(species))
        refuge <- refuge[sel_sp, , drop = FALSE]

        # Make data from from selected species
        plot_dat <- data.frame(w = rep(params@w, each = length(species)),
                               value = c(refuge),
                               Species = species)

        # Set up colors
        legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
        plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
        plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)
        linesize <- rep(0.8, length(legend_levels))
        names(linesize) <- names(params@linetype[legend_levels])

        p <- ggplot(plot_dat, aes(group = Species)) +
            facet_wrap(~ Species) +
            theme(strip.text.x = element_text(size = 6))
            #    strip.background = element_blank(),
            #   strip.text.x =element_blank())

        p + geom_line(aes(x = w, y = value,
                          colour = Legend, linetype = Legend,
                          linewidth = Legend)) +
            labs(colour = 'Functional Group', linetype = 'Functional Group',
                 linewidth = 'Functional Group') +
            scale_x_continuous(name = "Log Size [g]", trans = "log10",
                               breaks = c(10^-2, 10^0, 10^2, 10^4),
                               labels = c(-2, 0, 2, 4)) +
            scale_y_continuous(name = "Proportion Protected", limits = c(0, 1)) +
            scale_colour_manual(values = params@linecolour[legend_levels]) +
            scale_linetype_manual(values = params@linetype[legend_levels]) +
            scale_discrete_manual("linewidth", values = linesize)
    }
}


#' @rdname plotVulnerable
#' @export
plotlyRefuge <- function(object,
                             species = NULL,
                             ...) {

    argg <- as.list(environment())
    ggplotly(do.call("plotRefuge", argg),
             tooltip = c("Species", "w", "value"))
}


#' Plot the productivity of each species
#'
#' When called with a \linkS4class{MizerParams}
#' object the steady state productivity is plotted.
#'
#' @param param An object of class \linkS4class{MizerParams}
#' @param species The species to be selected. Optional. By default all
#'   species are selected. A vector of species names, or a numeric vector with
#'   the species indices, or a logical vector indicating for each species
#'   whether it is to be selected (TRUE) or not.
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept plots
#' @seealso [plotting_functions], [setRefuge()]
plotProductivity <- function(params,
                             species = NULL, time_range, highlight = NULL) {

    if (is(object, "MizerSim")) {

        stop('This functionality is not set up yet you dumbass
             .')
    } else if (is(object, "MizerParams")) {

        prod <- getProductivity(params)

        # selector for desired species
        sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                    error_on_empty = TRUE)
        species <- dimnames(params@initial_n)$sp[sel_sp]
        species <- gsub('inverts', NA, species)
        species <- species[!is.na(species)]
        sel_sp <- which(!is.na(species))
        prod <- prod[sel_sp, , drop = FALSE]

        # Make data from from selected species
        plot_dat <- data.frame(value = prod,Species = species)

        # Set up colors
        legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
        plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
        linesize <- rep(0.8, length(legend_levels))
        names(linesize) <- names(params@linetype[legend_levels])

        p <- ggplot(plot_dat, aes(group = Species))

        p + geom_bar(aes(x = Species, y = value,
                         colour = Legend, linetype = Legend,
                         linewidth = Legend)) +
            scale_y_continuous(name = "Productivity") +
            scale_colour_manual(values = params@linecolour[legend_levels])
    }
}

#' @rdname plotProductivity
#' @export
plotlyProductivity <- function(object,
                             species = NULL,
                             ...) {

    argg <- as.list(environment())
    ggplotly(do.call("plotProductivity", argg),
             tooltip = c("Species", "value"))
}

