#' Description of mizerReef plotting functions
#'
#' In addition to the plotting functions offered by the mizer package,
#' MizerReef provides new and extended plotting functions to visualize
#' both input parameters and the results of dynamic simulations for reef
#' ecosystem models. Several functions adapt or build on mizer's
#' plotting tools for reef-specific features.
#'
#' Available plotting functions:
#'
#' \strong{Input parameter plots}
#' \tabular{ll}{
#'   Plot \tab Description \cr
#'   [plotRefugeProfile()] \tab Plots the proportion of individuals (by length) that are protected from predators for each species (the refuge profile) \cr
#'   [plotDegradationScale()] \tab Plots a heatmap of the degradation scaling for refuge density in degradation simulations across bleaching years and size bins. \cr
#'   [plotDegScale()] \tab Plots a faceted heatmap comparing all three built-in degradation trajectories (rubble, algae, recovery) side by side. \cr
#'   [plotVulnerable()] \tab Plots vulnerability to predation by size and proportion for each species, either at steady state or for a chosen simulation time step. \cr
#'   [plotRefugeDensity()] \tab Plots the density of refuge by size at steady state and through time. \cr
#' }
#'
#' \strong{Result plots}
#' \tabular{ll}{
#'   Plot \tab Description \cr
#'   [mizer::plotBiomass()] \tab Plots the biomass of species and unstructured components through time (mizerReef extends this via a `getBiomass.mizerReefSim` method). \cr
#'   [mizer::plotSpectra2()] \tab Compare two size spectra (e.g., models or scenarios) in one plot. \cr
#'   [mizer::plotSpectraRelative()] \tab Plots relative difference between two spectra. \cr
#'   [plotSpectraChange()] \tab Plots the change (percent or relative proportion) between two spectra. \cr
#'   [plotRelativeContribution()] \tab Relative contribution of each species group to total abundance, biomass, and productivity. \cr
#'   [plotProductivity()] \tab Plots productivity for each species through time or at steady state. \cr
#'   [plot2Productivity()] \tab Productivity of two models or two size ranges in one plot. \cr
#'   [plotProductivityRelative()] \tab Plots the relative or percent change in productivity between two models or scenarios. \cr
#'   [plotTotalAbundance()] \tab Plots total abundance for each species at steady state. \cr
#'   [plotTotalBiomass()] \tab Plots total biomass for each species within a given size range at steady state. \cr
#'   [plot2TotalBiomass()] \tab Total biomass of two models or two size ranges in one plot. \cr
#'   [plotTotalBiomassRelative()] \tab Relative change in total biomass between two models or scenarios. \cr
#' }
#'
#' All plotting functions return a ggplot object, which can be further
#' customized using the ggplot2 grammar. Most functions accept either a
#' [MizerSim] or [MizerParams] object as input. Species colors and line
#' types are controlled by the `linecolour` and `linetype` slots in
#' [MizerParams], and can be changed by the user.
#'
#' For group or legend naming, you can add a column `group_names`
#' to your `species_params` data frame. This column should contain a
#' nicely formatted character string for each group or species, and will
#' be used in plot legends and facet labels for improved readability.
#'
#' Most plots allow selection of a subset of species via the `species`
#' argument. The order of species in legends matches the species parameter
#' data frame.
#'
#' @references
#' For the original mizer plotting functions and further details, see:
#' https://sizespectrum.org/mizer/reference/plotting_functions.html
#'
#' @section Usage:
#' Call these functions with a [MizerSim] or [MizerParams] object containing
#' simulation results. See individual function documentation for details and
#' examples.
#'
#' @author cmbeese
#' @name reef_plots
NULL

# Set global variables to avoid R CMD check notes
# Variables used - ggplot2 - known bug
# Hackiness to get past the 'no visible binding ... ' warning when running check
utils::globalVariables(c(
    "Species", "value", "Model", "Legend",
    "value.y", "value.x", "rel_diff", "l",
    "y_ticks", "highlight", "Metric",
    # variables used by ggplot for detritus and algae
    "Rate", "Source", "Consumer", "Year",
    "refuge_density", "size_bin", "scale_color_viridis",
    # variables used by ggplot in plotProductivity/plotRefugeDensity
    "Productivity", "time"
))

#' Plot the change between two spectra
#'
#' This plots the change between the steady state spectra of two mizer
#' objects. Let the spectra of the two objects be represented as
#' \eqn{N_1(w)} and \eqn{N_2(w)}. This function plots
#' \deqn{ \frac{N_2(w) - N_1(w)}{N_1(w)}}
#' expressed as a percentage (multiplied by 100) when `use_percent = TRUE`
#' (the default), or as the raw relative proportion when `use_percent =
#' FALSE`.
#'
#' For the difference calculated relative to the average of the two spectra,
#' \eqn{2 (N_2(w) - N_1(w)) / (N_2(w) + N_1(w))}, use mizer's own
#' [mizer::plotSpectraRelative()], which already dispatches correctly for
#' `mizerReef` objects.
#'
#' The individual spectra are calculated by the [mizer::plotSpectra()]
#' function which is passed all additional arguments you supply. So you can
#' for example determine a size range over which to average the simulation
#' results via the `time_range` argument. See [mizer::plotSpectra()] for more
#' options.
#'
#' @param object1 An object of class MizerSim or MizerParams
#' @param object2 An object of class MizerSim or MizerParams
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a
#'                  numeric vector with the species indices, or a logical
#'                  vector indicating for each species whether it is to be
#'                  selected (TRUE) or not.
#' @param power The abundance is plotted as the number density times the weight
#'   raised to this power. The default power = 1 gives the biomass density,
#'   whereas power = 2 gives the biomass density with respect to logarithmic
#'   size bins.
#' @param use_percent Logical. If TRUE (default), the change is expressed as
#'   a percentage (e.g. 50 for a 50% increase). If FALSE, the raw relative
#'   proportion is plotted instead (e.g. 0.5).
#' @param return_data Logical. If TRUE, returns the data frame underlying the
#'   plot instead of the plot itself. Default FALSE.
#' @param ... Parameters passed to `plotSpectra()`
#' @concept sumplots
#' @family plotting functions
#' @return A ggplot2 object, or a data frame if `return_data = TRUE`.
#' @export
plotSpectraChange <- function(object1, object2, species = NULL,
                              power, use_percent = TRUE,
                              return_data = FALSE, ...) {
    sf1 <- mizer::plotSpectra(object1,
        power = power,
        species = species,
        return_data = TRUE, ...
    )
    sf2 <- mizer::plotSpectra(object2,
        power = power,
        species = species,
        return_data = TRUE, ...
    )
    # plotSpectra() names the value column after the chosen power (e.g.
    # "Biomass density"); normalise it so the join below is stable.
    names(sf1)[2] <- "value"
    names(sf2)[2] <- "value"

    multiplier <- if (isTRUE(use_percent)) 100 else 1
    sf <- dplyr::left_join(sf1, sf2, by = c("w", "Legend")) |>
        dplyr::mutate(rel_diff = multiplier * (value.y - value.x) / value.x)
    yLabel <- if (isTRUE(use_percent)) "% Change in Biomass" else "Relative Change in Biomass"

    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }

    # group names -
    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }

    # Ensures that and species stay in the order we expect regardless of
    # naming
    sf$Legend <- factor(sf$Legend,
        levels = c(params@species_params$species, "Resource")
    )

    if (isTRUE(return_data)) {
        return(sf)
    }

    min_size <- min(sf$w)
    max_size <- max(sf$w)

    legend_levels <- intersect(
        names(params@linecolour),
        unique(sf$Legend)
    )
    linecolours <- params@linecolour[legend_levels]

    ggplot(sf, aes(x = w, y = rel_diff, colour = Legend)) +
        geom_line(linewidth = 0.95) +
        labs(x = "Weight [g]", y = yLabel) +
        scale_x_log10(limits = c(min_size, max_size)) +
        scale_color_manual(
            values = linecolours,
            labels = group_names[legend_levels]
        ) +
        geom_hline(
            yintercept = 0, linetype = 1,
            colour = "dark grey", linewidth = 1
        )
}

#' @importFrom plotly ggplotly
#' @rdname plotSpectraChange
#' @export
plotlySpectraChange <- function(object1, object2, ...) {
    ggplotly(plotSpectraChange(object1, object2, ...),
        tooltip = c("Legend", "w", "rel_diff")
    )
}


#' Plot the total productivity for each species Group
#'
#' When called with a [MizerParams] object the total steady
#' state productivity is plotted for each group. When called with a
#' [MizerSim] object the productivity of each species
#' through time is plotted.
#'
#' @inheritSection getProductivity Potential fisheries productivity
#'
#' @param object An object of class [MizerParams] or
#'                  [MizerSim]
#'
#' @param species The species to be selected. Optional. By default all species
#'   are selected. A vector of species names, or a numeric vector with the
#'   species indices, or a logical vector indicating for each species whether
#'   it is to be selected (TRUE) or not.
#'
#' @param start_time    The first time to be plotted. Default is the beginning
#'                      of the time series.
#'
#' @param end_time  The last time to be plotted. Default is the end of the time
#'                  series.
#'
#' @param facet   A boolean value indicating whether to facet the result plot by
#'                species group. Defaults to TRUE.
#'
#' @param total   A boolean value that determines whether the total productivity
#'                from all species is plotted as well. Default is TRUE.
#'
#' @param include_inverts   A boolean value that determines whether the
#'                      "inverts" species group is included. Default is
#'                      FALSE, since invertebrate productivity is typically
#'                      not relevant to fishing yield. Only takes effect
#'                      when `species` is not explicitly provided.
#'
#' @param return_data   A boolean value that determines whether the formatted
#'                      data used for the plot is returned instead of the
#'                      plot itself. Default value is FALSE.
#'
#' @inheritParams getProductivity
#'
#' @param ... unused
#'
#' @return  A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'          frame with the the productivity for each Species Group
#'          is returned.
#'
#' @import ggplot2
#' @export
#'
#' @concept sumplots
#' @family plotting functions
#'
#' @seealso [getEGrowthTime()],[getProductivity()],
#'          [plotBiomass()], [plot2TotalBiomass()],
#'          [plotTotalBiomassRelative()], [plotProductivity()],
#'          [plot2Productivity()], [plotProductivityRelative()]
plotProductivity <- function(object,
                             start_time = NULL, end_time = NULL,
                             facet = TRUE, species = NULL, total = TRUE,
                             min_fishing_l = NULL, max_fishing_l = NULL,
                             include_repro = FALSE, include_inverts = FALSE,
                             return_data = FALSE, ...) {
    if (is(object, "MizerSim")) {
        # sim values
        sim <- object
        params <- sim@params
        assert_that(
            is(sim, "MizerSim"),
            is.flag(return_data)
        )

        if (is.null(species)) {
            # is.null(), not missing(): wrapper functions like
            # plot2Productivity()/plotProductivityRelative() always forward
            # species = species (their own default is also NULL), which
            # makes missing() FALSE even when the top-level caller never
            # specified a species -- silently breaking both the "all
            # species" default and the include_inverts default below
            # whenever plotProductivity() is called via one of those
            # wrappers rather than directly.
            species <- params@species_params$species
            if (!isTRUE(include_inverts)) {
                species <- setdiff(species, "inverts")
            }
        }

        if (missing(start_time)) {
            start_time <- as.numeric(dimnames(sim@n)[[1]][1])
        }
        if (missing(end_time)) {
            end_time <- as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]])
        }
        if (start_time >= end_time) {
            stop("start_time must be less than end_time")
        }

        time_range <- start_time:end_time

        p <- getProductivity(sim,
            time_range = time_range,
            min_fishing_l = min_fishing_l,
            max_fishing_l = max_fishing_l,
            include_repro = include_repro, ...
        )

        # Include total
        if (total) {
            p <- cbind(p, Total = rowSums(p))
        }

        p <- reshape2::melt(p)

        # Implement ylim and a minimal cutoff and bring columns in
        # desired order
        names(p) <- c("Year", "Species", "Productivity")

        # Select species
        plot_dat <- p[p$Species %in% c("Total", species), ]
        plot_dat$Legend <- plot_dat$Species
        legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
        plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)

        if (return_data) {
            return(plot_dat)
        }

        p <- ggplot(plot_dat, aes(
            x = Year, y = Productivity,
            group = Legend, colour = Legend,
            linetype = Legend
        ))

        if (facet == TRUE) {
            p + facet_wrap(~Legend, scales = "free_y")
        }

        p <- p +
            geom_line(linewidth = 0.8) +
            scale_colour_manual(values = params@linecolour[legend_levels]) +
            scale_linetype_manual(values = params@linetype[legend_levels]) +
            labs(
                colour = "Species Group",
                linetype = "Species Group",
                x = "Year"
            )
    } else if (is(object, "MizerParams")) {
        # params
        params <- object
        assert_that(
            is(params, "MizerParams"),
            is.flag(return_data)
        )

        ### group names
        if (is.null(params@species_params$group_names)) {
            group_names <- params@species_params$species
            names(group_names) <- params@species_params$species
        } else {
            group_names <- params@species_params$group_names
            names(group_names) <- params@species_params$species
        }

        ### get productivity
        prod <- getProductivity(params,
            min_fishing_l = min_fishing_l,
            max_fishing_l = max_fishing_l,
            include_repro = include_repro
        )

        ### species selector
        sel_sp <- valid_species_arg(params, species,
            return.logical = TRUE,
            error_on_empty = TRUE
        )
        # sel_sp stays a logical vector over the *original* species order
        # throughout, so it can subset both `species` and `prod`
        # consistently -- recomputing it from an already-filtered vector
        # (as a previous version of this code did) silently misaligns
        # species labels with the wrong data once the excluded species
        # isn't last in the species order.
        if (is.null(species) && !isTRUE(include_inverts)) {
            # is.null(), not missing() -- see the matching comment in the
            # MizerSim branch above.
            sel_sp <- sel_sp & (dimnames(params@initial_n)$sp != "inverts")
        }
        species <- dimnames(params@initial_n)$sp[sel_sp]
        prod <- prod[sel_sp, drop = FALSE]

        ### data frame from selected species
        plot_dat <- data.frame(value = prod, Species = species)

        ### colors
        legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
        plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)

        ### return data if requested -
        if (return_data) {
            return(plot_dat)
        }

        # plot -
        p <- ggplot(plot_dat, aes(
            x = Species, y = value,
            group = Legend, fill = Legend
        ))

        p + geom_bar(stat = "identity", position = "dodge") +
            scale_y_continuous(name = expression(Productivity ~ "(" * g / m^2 / year * ")")) +
            scale_fill_manual(
                values = params@linecolour[legend_levels],
                labels = group_names[legend_levels]
            ) +
            labs(fill = "Species Group", x = "Species Group")
    } else {
        stop("Object should be a MizerParams or MizerSim object.")
    }
}


#' @importFrom plotly ggplotly
#' @rdname plotProductivity
#' @export
plotlyProductivity <- function(object,
                               species = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plot2Productivity", argg),
        tooltip = c("Species", "value")
    )
}


#' Plot the fisheries productivity of two models or two different size ranges
#'  in the same plot
#'
#' When called with a [MizerParams]
#' object the steady state productivities are plotted.
#'
#' @param object1 First MizerParams or MizerSim object.
#'
#' @param object2 Second MizerParams or MizerSim object.
#'
#' @param species   The species to be selected. Optional. By default all target
#'                  species are selected. A vector of species names, or a numeric
#'                  vector with the species indices, or a logical vector
#'                  indicating for each group whether it is to be selected
#'                  (TRUE) or not.
#'
#' @param name1 An optional string with the name for the first model, to be
#'              used in the legend. Set to "First" by default.
#'
#' @param name2 An optional string with the name for the second model, to be
#'              used in the legend. Set to "Second" by default.
#'
#' @param min_fishing_l1    Optional.  The minimum length (cm) of fished
#'                          individuals for model 2. Defaults to 7cm.
#'                          A parameter passed to [getProductivity()].
#'
#' @param max_fishing_l1    Optional.  The maximum length (cm) of fished
#'                          individuals for model 1. Defaults to max length.
#'                          A parameter passed to [getProductivity()].
#'
#' @param min_fishing_l2    Optional.  The minimum length (cm) of fished
#'                          individuals for model 2. Defaults to 7cm.
#'                          A parameter passed to [getProductivity()].
#'
#' @param max_fishing_l2    Optional.  The maximum length (cm) of fished
#'                          individuals for model 1. Defaults to max length.
#'                          A parameter passed to [getProductivity()].
#'
#' @param stack     A boolean value that determines whether bars are separated
#'                  by species. Defaults to FALSE. If true, returns a stacked
#'                  barplot with the total biomass for each group instead of
#'                  individual bars for each group. Useful for comparison
#'                  between steady states.
#'
#' @param return_data   A boolean value that determines whether the formatted
#'                      data used for the plot is returned instead of the plot
#'                      itself. Default value is FALSE.
#'
#' @inheritDotParams plotProductivity
#'
#' @return  A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'          frame with the the productivity for each Species Group
#'          by model is returned.
#'
#' @import ggplot2
#' @export
#'
#' @concept sumplots
#' @family plotting functions
#'
#' @seealso [plotBiomass()], [plot2TotalBiomass()], [plotTotalBiomassRelative()],
#'          [plotProductivity()], [plot2Productivity()], [plotProductivityRelative()]
plot2Productivity <- function(object1, object2, species = NULL,
                              name1 = "First", name2 = "Second",
                              min_fishing_l1 = NULL, max_fishing_l1 = NULL,
                              min_fishing_l2 = NULL, max_fishing_l2 = NULL,
                              stack = FALSE,
                              return_data = FALSE, ...) {
    # get data frames with plotProductivity -
    sf1 <- plotProductivity(object1,
        species = species,
        drop = TRUE,
        min_fishing_l = min_fishing_l1,
        max_fishing_l = max_fishing_l1,
        return_data = TRUE, ...
    )
    sf1$Model <- name1
    sf2 <- plotProductivity(object2,
        species = species,
        drop = TRUE,
        min_fishing_l = min_fishing_l2,
        max_fishing_l = max_fishing_l2,
        return_data = TRUE, ...
    )
    sf2$Model <- name2

    sf <- rbind(sf1, sf2)

    # Make sure order of models isnt changed by names
    sf$Model <- factor(sf$Model, levels = c(name1, name2))

    # if sim, get params -
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }

    # group names -
    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }

    # plot -
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    sf$Legend <- factor(sf$Species, levels = legend_levels)

    # Return data frame if requested
    if (return_data == TRUE) {
        return(sf)
    }

    if (stack == FALSE) {
        p <- ggplot(sf, aes(
            x = Species, y = value,
            group = Model, alpha = Model, fill = Legend
        ))

        p + geom_bar(stat = "identity", position = "dodge", color = "black") +
            scale_y_continuous(name = expression(Productivity ~ "(" * g / m^2 / year * ")")) +
            scale_fill_manual(
                values = params@linecolour[legend_levels],
                labels = group_names[legend_levels]
            ) +
            scale_alpha_manual(
                values = c(0.5, 1),
                labels = c(name1, name2)
            ) +
            labs(fill = "Species Group", x = "Species Group")
    } else if (stack == TRUE) {
        p <- ggplot(sf, aes(
            x = Model, y = value,
            alpha = Model, fill = Legend
        ))

        p + geom_bar(stat = "identity", position = "stack", color = "black") +
            scale_y_continuous(name = expression(Productivity ~ "(" * g / m^2 / year * ")")) +
            scale_fill_manual(
                values = params@linecolour[legend_levels],
                labels = group_names[legend_levels]
            ) +
            scale_alpha_manual(
                values = c(0.5, 1),
                labels = c(name1, name2)
            ) +
            labs(fill = "Species Group", x = "Model")
    }
}

#' @importFrom plotly ggplotly
#' @rdname plot2Productivity
#' @export
plotly2Productivity <- function(object1, object2, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plot2Productivity", argg),
        tooltip = c("Species", "value")
    )
}

#' Plot the relative difference between the potential fisheries productivity
#' rates of two models or two different size ranges in the same plot
#'
#' This function creates a barplot with the percent change in potential
#' fisheries productivity between either: (1) two different mizerParams objects
#' (2) two different size ranges. If comparing two mizerParams objects, they
#' must have the same species.
#'
#' This function is usually used in conjunction with
#' [plotTotalBiomassRelative()] to check for decoupling between biomass and
#' productivity.
#'
#' The individual productivity rates are calculated by the
#' [plotProductivity()] function which is passed all additional arguments you
#' supply. See [plotProductivity()] for more details.
#'
#' To compare between different size ranges, use the `min_fishing_l1`
#' and `max_fishing_l1` arguments for the first size range and  the
#' `min_fishing_l2`and `max_fishing_l2` arguments for the second.
#'
#' @param object1 First MizerParams or MizerSim object.
#'
#' @param object2 Second MizerParams or MizerSim object.
#'
#' @param species   The species to be selected. Optional. By default all target
#'                  species are selected. A vector of species names, or a numeric
#'                  vector with the species indices, or a logical vector
#'                  indicating for each species whether it is to be selected
#'                  (TRUE) or not.
#'
#' @param diff_method   The method to calculate the relative change between
#'                      models. If `percent_change`, the percent change is
#'                      calculated relative to the value from object 1 with
#'                      formula 100*(new-old)/old. If `rel_diff` the relative
#'                      difference is returned given by (new - old)/(old + new).
#'
#' @param min_fishing_l1    Optional.  The minimum length (cm) of fished
#'                          individuals for model 2. Defaults to 7cm.
#'                          A parameter passed to [getProductivity()].
#'
#' @param max_fishing_l1    Optional.  The maximum length (cm) of fished
#'                          individuals for model 1. Defaults to max length.
#'                          A parameter passed to [getProductivity()].
#'
#' @param min_fishing_l2    Optional.  The minimum length (cm) of fished
#'                          individuals for model 2. Defaults to 7cm.
#'                          A parameter passed to [getProductivity()].
#'
#' @param max_fishing_l2    Optional.  The maximum length (cm) of fished
#'                          individuals for model 1. Defaults to max length.
#'                          A parameter passed to [getProductivity()].
#'
#' @param return_data   A boolean value that determines whether the formatted
#'                      data used for the plot is returned instead of the plot
#'                      itself. Default value is FALSE.
#'
#' @inheritDotParams plotProductivity
#'
#' @return  A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'          frame with the the productivity for each species
#'          by model is returned as well as another column called `rel_diff`
#'          that gives the relative difference between the two values.
#'
#' @import ggplot2
#' @export
#'
#' @concept sumplots
#' @family plotting functions
#'
#' @seealso [plotBiomass()], [plot2TotalBiomass()], [plotTotalBiomassRelative()],
#'          [plotProductivity()], [plot2Productivity()], [plotProductivityRelative()]
plotProductivityRelative <- function(object1, object2, diff_method,
                                     species = NULL,
                                     min_fishing_l1 = NULL,
                                     max_fishing_l1 = NULL,
                                     min_fishing_l2 = NULL,
                                     max_fishing_l2 = NULL,
                                     return_data = FALSE, ...) {
    # get data frames with plotProductivity -
    sf1 <- plotProductivity(object1,
        species = species,
        min_fishing_l = min_fishing_l1,
        max_fishing_l = max_fishing_l1,
        return_data = TRUE, ...
    )
    sf2 <- plotProductivity(object2,
        species = species,
        min_fishing_l = min_fishing_l2,
        max_fishing_l = max_fishing_l2,
        return_data = TRUE, ...
    )

    # Calculate relative difference
    if (diff_method == "percent_change") {
        sf <- dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |>
            dplyr::mutate(rel_diff = 100 * (value.y - value.x) / value.x)
        yLabel <- "% Change in Productivity"
    } else if (diff_method == "rel_diff") {
        sf <- dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |>
            dplyr::mutate(rel_diff = ((value.y - value.x) / (value.x + value.y)))
        yLabel <- "Relative Difference in Productivity"
    } else {
        stop("diff_method should be either 'percent_change' or 'rel_diff'.")
    }

    # return data if requested -
    if (return_data == TRUE) {
        return(sf)
    }

    # save 1 set for species names -
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }

    # species names -
    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }

    # plot --
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    sf$Legend <- factor(sf$Species, levels = legend_levels)

    p <- ggplot(sf, aes(
        x = Species, y = rel_diff,
        fill = Legend
    ))

    p + geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_y_continuous(name = yLabel) +
        scale_fill_manual(
            values = params@linecolour[legend_levels],
            labels = group_names[legend_levels]
        ) +
        labs(fill = "Species", x = "Species") +
        geom_hline(
            yintercept = 0, linetype = 1,
            colour = "dark grey", linewidth = 0.9
        )
}

#' @importFrom plotly ggplotly
#' @rdname plotProductivityRelative
#' @export
plotlyProductivityRelative <- function(object1, object2,
                                       diff_method, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotProductivityRelative", argg),
        tooltip = c("Species", "value")
    )
}

#' Plot the total abundance for each species in a size range
#'
#' This functions creates a barplot with the abundance of each species
#' within a given size range.
#'
#' @param object An object of class [MizerParams] or [MizerSim].
#'               If a [MizerSim] object is provided, the abundance at the
#'              last time step is used. 
#'
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a
#'                  numeric vector with the species indices, or a logical
#'                  vector indicating for each species whether it is to be
#'                  selected (TRUE) or not.
#'
#' @param min_fishing_l The minimum length (cm) for abundance estimates.
#'                      Defaults to smallest size.
#'
#' @param max_fishing_l The maximum length (cm) of for abundance estimates.
#'                      Defaults to max length.
#'
#' @param return_data   A boolean value that determines whether the formatted
#'                      data used for the plot is returned instead of the plot
#'                      itself. Default value is FALSE.
#'
#' @param stack     A boolean value that determines whether bars are separated
#'                  by species. Defaults to FALSE. If true, returns a stacked
#'                  barplot with the total biomass for each species instead of
#'                  individual bars for each species. Useful for comparison
#'                  between steady states.
#'
#' @param ... unused
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @export
#'
#' @family plotting functions
#' @concept sumplots
#' @seealso [plotBiomass()], [plot2TotalBiomass()],
#'          [plotTotalBiomassRelative()],
#'          [plotProductivity()],
#'          [plot2Productivity()], [plotProductivityRelative()]
plotTotalAbundance <- function(object, stack = FALSE,
                               species = NULL,
                               min_fishing_l = NULL,
                               max_fishing_l = NULL,
                               return_data = FALSE, ...) {
    # object checks -
    if (is(object, "MizerSim")) {
        ## sim values -
        # get total abundance at last timestep
        params <- object@params

        abd <- mizer::getN(object,
            min_l = min_fishing_l,
            max_l = max_fishing_l
        )

        # Use the last *row position*, not the numeric time value as an
        # index -- row 1 is time 0, so treating the time value itself as a
        # positional index is off by one (e.g. abd[2, ] would grab the row
        # for time = 1, not time = 2). Matches the correct approach already
        # used in plotProductivity()'s own MizerSim branch.
        abd <- abd[dim(abd)[1], , drop = TRUE]

    } else {
        # params -
        params <- object
        assert_that(
            is(params, "MizerParams"),
            is.flag(return_data)
        )

        abd <- mizer::getN(params,
            min_l = min_fishing_l,
            max_l = max_fishing_l
        )
    }

    # create plot_dat -
    ## group names -
    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }

    ## species selector -
    sel_sp <- mizer::valid_species_arg(params, species,
        return.logical = TRUE,
        error_on_empty = TRUE
    )
    species <- dimnames(params@initial_n)$sp[sel_sp]
    abd <- abd[sel_sp, drop = FALSE]

    ## data frame from selected species -
    plot_dat <- data.frame(value = abd, Species = species)

    ## colors -
    legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)

    ## return data if requested -
    if (return_data) {
        return(plot_dat)
    }

    if (stack == FALSE) {
        p <- ggplot(plot_dat, aes(x = Species, y = value, fill = Legend))

        p + geom_bar(stat = "identity", position = "dodge") +
            scale_y_continuous(name = expression("Total Abundance (no./m^2)")) +
            scale_fill_manual(
                values = params@linecolour[legend_levels],
                labels = group_names[legend_levels]
            )
    } else if (stack == TRUE) {
        p <- ggplot(plot_dat, aes(x = Species, y = value, fill = Legend))

        p + geom_bar(stat = "identity", position = "stack") +
            scale_y_continuous(name = expression("Total Abundance (no./m^2)")) +
            scale_fill_manual(
                values = params@linecolour[legend_levels],
                labels = group_names[legend_levels]
            )
    }
}

#' @importFrom plotly ggplotly
#' @rdname plotTotalAbundance
#' @export
plotlyTotalAbundance <- function(object, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotTotalAbundance", argg),
        tooltip = c("Species", "value")
    )
}


#' Plot the total biomass for each species in a size range
#'
#' This functions creates a barplot with the biomass of each specie
#' within a size range. Usually used in conjunction with [plotProductivity()]
#' to check for decoupling.
#'
#' @param object An object of class [MizerParams] or [MizerSim].
#'               If a [MizerSim] object is provided, the biomass at the
#'               last time step is used. 
#'
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a
#'                  numeric vector with the species indices, or a logical
#'                  vector indicating for each species whether it is to be
#'                  selected (TRUE) or not.
#'
#' @param min_fishing_l The minimum length (cm) for biomass estimates.
#'                      Defaults to smallest size.
#'
#' @param max_fishing_l The maximum length (cm) of for biomass estimates.
#'                      Defaults to max length.
#'
#' @param return_data   A boolean value that determines whether the formatted
#'                      data used for the plot is returned instead of the plot
#'                      itself. Default value is FALSE.
#'
#' @param ... unused
#'
#' @return A ggplot2 object
#'
#' @import ggplot2
#' @export
#'
#' @family plotting functions
#' @concept sumplots
#' @seealso [plotBiomass()], [plot2TotalBiomass()],
#'          [plotTotalBiomassRelative()],
#'          [plotProductivity()], [plot2Productivity()],
#'          [plotProductivityRelative()]
plotTotalBiomass <- function(object,
                             species = NULL,
                             min_fishing_l = NULL,
                             max_fishing_l = NULL,
                             return_data = FALSE, ...) {
    # object checks -
    if (is(object, "MizerSim")) {
        ## sim values -
        # get total biomass at last timestep
        params <- object@params

        biom <- mizer::getBiomass(object,
            min_l = min_fishing_l,
            max_l = max_fishing_l
        )

        # Use the last *row position*, not the numeric time value as an
        # index -- see the matching comment in plotTotalAbundance().
        biom <- biom[dim(biom)[1], , drop = TRUE]

    } else {
        # params -
        params <- object
        assert_that(
            is(params, "MizerParams"),
            is.flag(return_data)
        )

        biom <- mizer::getBiomass(params,
            min_l = min_fishing_l,
            max_l = max_fishing_l
        )
    }

    # create plot_dat -

    ## group names -
    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }

    ## species selector -
    # sel_sp stays a logical vector over the *original* species order so it
    # subsets `species` and `biom` consistently -- a previous version of
    # this code recomputed sel_sp <- which(!is.na(species)) on the already-
    # sliced `species` vector, which is always a no-op is.na() check with
    # nothing upstream ever introducing an NA, so it silently collapsed to
    # 1:length(species) and misaligned biom values with the wrong species
    # labels whenever `species` selected a non-contiguous subset (e.g.
    # species = c("predators", "inverts"), skipping "herbivores" in the
    # middle of the original species order).
    sel_sp <- mizer::valid_species_arg(params, species,
        return.logical = TRUE,
        error_on_empty = TRUE
    )
    species <- dimnames(params@initial_n)$sp[sel_sp]
    # Index by name, not by the sel_sp logical vector directly: for a
    # MizerSim, biom comes from getBiomass.mizerReefSim(), which also
    # includes algae/detritus columns alongside the species -- indexing
    # that longer, differently-shaped vector with a species-only logical
    # vector recycles it incorrectly instead of erroring.
    biom <- biom[species, drop = FALSE]

    ## data frame from selected species -
    plot_dat <- data.frame(value = biom, Species = species)

    ## colors -
    legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)

    ## return data if requested -
    if (return_data) {
        return(plot_dat)
    }

    # plot -
    p <- ggplot(plot_dat, aes(
        x = Species, y = value,
        group = Legend, fill = Legend
    ))

    p + geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(name = expression("Total Biomass" ~ "(" * g / m^2 * ")")) +
        scale_fill_manual(
            values = params@linecolour[legend_levels],
            labels = group_names[legend_levels]
        ) +
        labs(fill = "Species", x = "Species")
}

#' @importFrom plotly ggplotly
#' @rdname plotTotalBiomass
#' @export
plotlyTotalBiomass <- function(object, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotTotalBiomass", argg),
        tooltip = c("Species", "value")
    )
}


#' Plot the total biomass of two models or of two different size ranges in
#' the same plot
#'
#' When called with a [MizerParams]
#' object the steady state biomasses are plotted.
#'
#' @inheritParams plot2Productivity
#'
#' @inheritDotParams plotTotalBiomass
#'
#' @return  A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'          frame with the the total steady state biomass for each species
#'          by model is returned as well as another column called
#'          `rel_diff`that gives the relative difference between the two
#'          values.
#'
#' @import ggplot2
#' @export
#'
#' @concept sumplots
#' @family plotting functions
#'
#' @seealso [plotBiomass()], [plot2TotalBiomass()], [plotTotalBiomassRelative()],
#'          [plotProductivity()], [plot2Productivity()], [plotProductivityRelative()]
plot2TotalBiomass <- function(object1, object2,
                              species = NULL,
                              name1 = "First", name2 = "Second",
                              min_fishing_l1 = NULL, max_fishing_l1 = NULL,
                              min_fishing_l2 = NULL, max_fishing_l2 = NULL,
                              stack = FALSE,
                              return_data = FALSE, ...) {
    # get data frames with plotTotalBiomass -
    sf1 <- plotTotalBiomass(object1,
        species = species,
        min_fishing_l = min_fishing_l1,
        max_fishing_l = max_fishing_l1,
        return_data = TRUE, ...
    )
    sf1$Model <- name1
    sf2 <- plotTotalBiomass(object2,
        species = species,
        min_fishing_l = min_fishing_l2,
        max_fishing_l = max_fishing_l2,
        return_data = TRUE, ...
    )
    sf2$Model <- name2

    sf <- rbind(sf1, sf2)

    # Make sure model names dont change order
    sf$Model <- factor(sf$Model, levels = c(name1, name2))

    # if sim, get params -
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }

    # group names -
    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }

    # plot -
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    sf$Legend <- factor(sf$Species, levels = legend_levels)

    # Return data frame if requested
    if (return_data == TRUE) {
        return(sf)
    }

    if (stack == FALSE) {
        p <- ggplot(sf, aes(
            x = Species, y = value,
            group = Model, alpha = Model, fill = Legend
        ))

        p + geom_bar(stat = "identity", position = "dodge", color = "black") +
            scale_y_continuous(name = expression("Total Biomass" ~ "(" * g / m^2 * ")")) +
            scale_fill_manual(
                values = params@linecolour[legend_levels],
                labels = group_names[legend_levels]
            ) +
            scale_alpha_manual(
                values = c(0.5, 1),
                labels = c(name1, name2)
            ) +
            labs(fill = "Species", x = "Species")
    } else if (stack == TRUE) {
        p <- ggplot(sf, aes(
            x = Model, y = value,
            alpha = Model, fill = Legend
        ))

        p + geom_bar(stat = "identity", position = "stack", color = "black") +
            scale_y_continuous(name = expression("Total Biomass" ~ "(" * g / m^2 * ")")) +
            scale_fill_manual(
                values = params@linecolour[legend_levels],
                labels = group_names[legend_levels]
            ) +
            scale_alpha_manual(
                values = c(0.5, 1),
                labels = c(name1, name2)
            ) +
            labs(fill = "Species", x = "Model")
    }
}

#' @importFrom plotly ggplotly
#' @rdname plot2TotalBiomass
#' @export
plotly2TotalBiomass <- function(object1, object2,
                                species = NULL, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plot2TotalBiomass", argg),
        tooltip = c("Species", "value")
    )
}

#' Plot the relative difference in between the total biomasses of each
#' each species within a size range at steady state
#'
#' This functions creates a barplot with the relative change in biomass of
#' each species within a size range between either (1) two different
#' mizerParams objects (two models) or (2) two different size ranges.
#'
#' This function is usually used in conjunction with
#' [plotProductivityRelative()] to check for decoupling between biomass and
#' productivity.
#'
#' The individual productivity rates are calculated by the
#' [plotTotalBiomass()] function which is passed all additional arguments you
#' supply. See [plotTotalBiomass()] for more details.
#'
#' To compare between different size ranges, use the `min_fishing_l1`
#' and `max_fishing_l1` arguments for the first size range and  the
#' `min_fishing_l2`and `max_fishing_l2` arguments for the second.
#'
#' @inheritParams plotProductivityRelative
#'
#' @inheritDotParams plotTotalBiomass
#'
#' @return  A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'          frame with the the total steady state biomass for each functional
#'          group by model is returned as well as another column called
#'          `rel_diff`that gives the relative difference between the two
#'          values.
#'
#' @import ggplot2
#' @export
#'
#' @concept sumplots
#' @family plotting functions
#'
#' @seealso [plotBiomass()], [plot2TotalBiomass()], [plotTotalBiomassRelative()],
#'          [plotProductivity()], [plot2Productivity()], [plotProductivityRelative()]
plotTotalBiomassRelative <- function(object1, object2,
                                     diff_method,
                                     species = NULL,
                                     min_fishing_l1 = NULL,
                                     max_fishing_l1 = NULL,
                                     min_fishing_l2 = NULL,
                                     max_fishing_l2 = NULL,
                                     return_data = FALSE, ...) {
    # get data frames with plotTotalBiomass -
    sf1 <- plotTotalBiomass(object1,
        species = species,
        min_fishing_l = min_fishing_l1,
        max_fishing_l = max_fishing_l1,
        return_data = TRUE, ...
    )
    sf2 <- plotTotalBiomass(object2,
        species = species,
        min_fishing_l = min_fishing_l2,
        max_fishing_l = max_fishing_l2,
        return_data = TRUE, ...
    )

    # Calculate relative difference
    if (diff_method == "percent_change") {
        sf <- dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |>
            dplyr::mutate(rel_diff = 100 * (value.y - value.x) / value.x)

        yLabel <- "% Change in Total Biomass"
    } else if (diff_method == "rel_diff") {
        sf <- dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |>
            dplyr::mutate(rel_diff = ((value.y - value.x) / (value.x + value.y)))

        yLabel <- "Relative Difference in Total Biomass"
    } else {
        stop("diff_method should be either 'percent_change' or 'rel_diff'.")
    }

    # Return data frame if requested
    if (return_data == TRUE) {
        return(sf)
    }

    # if sim, get params -
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }

    # group names -
    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }

    # plot -
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    sf$Legend <- factor(sf$Species, levels = legend_levels)

    p <- ggplot(sf, aes(x = Species, y = rel_diff, fill = Legend))

    p + geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_y_continuous(name = yLabel) +
        scale_fill_manual(
            values = params@linecolour[legend_levels],
            labels = group_names[legend_levels]
        ) +
        labs(fill = "Species Group", x = "Species Group") +
        geom_hline(
            yintercept = 0, linetype = 1,
            colour = "dark grey", linewidth = 0.9
        )
}

#' @importFrom plotly ggplotly
#' @rdname plotTotalBiomassRelative
#' @export
plotlyTotalBiomassRelative <- function(object1, object2,
                                       diff_method, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("TotalBiomassRelative", argg),
        tooltip = c("Species", "value")
    )
}


#' Plot the relative contribution of each species group to total abundance,
#' total biomass, and total productivity
#'
#' The group abundances, biomasses, productivities are calculated by the
#' [plotTotalAbundance()], [plotTotalBiomass()], and [plotProductivity()]
#'  functions. These are passed all additional arguments you supply. See
#'  [plotTotalAbundance()], [plotTotalBiomass()] and [plotProductivity()]
#'  for more details.
#'
#' @param object An object of class [MizerParams]
#'
#' @param min_size  parameters be passed to [plotTotalAbundance()] and
#'                  [plotTotalBiomass()]. The minimum length (cm) of
#'                  individuals for biomass estimates. Defaults to
#'                  smallest size in the model.
#'
#' @param min_fishing_l parameters be passed to [getProductivity()]. The
#'                      minimum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to 7 cm.
#'
#' @param return_data   A boolean value that determines whether the formatted
#'                      data used for the plot is returned instead of the plot
#'                      itself. Default value is FALSE.
#'
#' @inheritDotParams plotTotalBiomass
#' @inheritDotParams plotTotalAbundance
#' @inheritDotParams plotProductivity
#'
#' @param include_inverts   A boolean value that determines whether the
#'                      "inverts" species group is included. Default is
#'                      FALSE, since invertebrate productivity is typically
#'                      not relevant to fishing yield.
#'
#' @import ggplot2
#' @export
#'
#' @concept sumplots
#' @family plotting functions
#' @seealso [plotTotalAbundance()], [plotTotalBiomass()], [plotProductivity()]
plotRelativeContribution <- function(object,
                                     min_size = NULL,
                                     min_fishing_l = NULL,
                                     include_inverts = FALSE,
                                     return_data = FALSE, ...) {
    abd <- plotTotalAbundance(object,
        min_fishing_l = min_size,
        return_data = TRUE, ...
    )
    abd$Metric <- "Abundance"

    biom <- plotTotalBiomass(object,
        min_fishing_l = min_size,
        return_data = TRUE, ...
    )
    biom$Metric <- "Biomass"

    prod <- plotProductivity(object,
        min_fishing_l = min_fishing_l,
        include_inverts = include_inverts,
        return_data = TRUE, ...
    )
    prod$Metric <- "Productivity"

    if (is(object, "MizerSim")) {
        params <- object@params
    } else {
        params <- object
        assert_that(is(params, "MizerParams"))
    }

    # group names -
    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }

    # Remove invertebrates (prod is already filtered via plotProductivity()'s
    # own include_inverts argument above; abd/biom need it applied here since
    # plotTotalAbundance()/plotTotalBiomass() don't filter species by name)
    if (!isTRUE(include_inverts)) {
        abd <- subset(abd, Species != "inverts")
        biom <- subset(biom, Species != "inverts")
    }

    # Relative Contribution
    abd <- dplyr::mutate(abd, rel = value / sum(value))
    biom <- dplyr::mutate(biom, rel = value / sum(value))
    prod <- dplyr::mutate(prod, rel = value / sum(value))

    rel <- rbind(abd, biom, prod)

    # Legend
    legend_levels <- intersect(names(params@linecolour), unique(rel$Legend))
    rel$Legend <- factor(rel$Species, levels = legend_levels)
    rel$Metric <- factor(rel$Metric, levels = c(
        "Abundance",
        "Biomass",
        "Productivity"
    ))

    # Return data if requested
    if (return_data == TRUE) {
        return(rel)
    }

    # Plot
    p <- ggplot(rel, aes(x = rel, y = Metric, fill = Legend))

    p + geom_bar(stat = "identity", position = "fill") +
        scale_fill_manual(
            values = params@linecolour[legend_levels],
            labels = group_names[legend_levels]
        ) +
        labs(x = "Relative Contribution", y = "Metric", fill = "")
}

#' @importFrom plotly ggplotly
#' @rdname plotRelativeContribution
#' @export
plotlyRelativeContribution <- function(object, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotRelativeContribution", argg),
        tooltip = c("Species", "rel", "Metric")
    )
}

#' Plot refuge density through time
#'
#' Plots the refuge density at each size bin through time as colored lines,
#' showing how refuge degrades over the course of a simulation. Uses viridis
#' color scale to represent temporal succession.
#'
#' @param sim A [MizerSim] object
#' @param return_data A boolean value that determines whether the formatted
#'                    data used for the plot is returned instead of the plot
#'                    itself. Default value is FALSE.
#' @param ... Unused
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'         frame with refuge density at each size and time
#'
#' @import ggplot2
#' @export
#' @family plotting functions
#' @concept refugePlots
#' @seealso [plotting_functions], [setRefuge()], [reefDegrade()]
plotRefugeDensity <- function(sim, return_data = FALSE, ...) {
    # Input validation
    assertthat::assert_that(
        inherits(sim, "MizerSim"),
        is.logical(return_data),
        length(return_data) == 1
    )

    # Use getDegrade to get refuge density array (time x size_bin)
    degrade <- getDegrade(sim)

    plot_dat <- as.data.frame(as.table(degrade))
    colnames(plot_dat) <- c("time", "size_bin", "refuge_density")
    plot_dat$time <- as.numeric(as.character(plot_dat$time))
    plot_dat$size_bin <- as.numeric(as.character(plot_dat$size_bin))

    # Return data if requested
    if (isTRUE(return_data)) {
        return(plot_dat)
    }

    # Create plot
    p <- ggplot(plot_dat, aes(
        x = time, y = refuge_density,
        group = factor(size_bin),
        color = factor(size_bin)
    )) +
        geom_line(linewidth = 0.8) +
        scale_color_brewer(
            name = "Size Bin (cm)",
            palette = "Spectral"
        ) +
        labs(
            x = "Time (years)",
            y = "Refuge Density",
            title = "Refuge Density Through Time"
        ) +
        theme_minimal() +
        theme(legend.position = "right")

    return(p)
}

#' @importFrom plotly ggplotly
#' @rdname plotRefugeDensity
#' @export
plotlyRefugeDensity <- function(sim, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotRefugeDensity", argg),
        tooltip = c("time", "size_bin", "refuge_density")
    )
}
