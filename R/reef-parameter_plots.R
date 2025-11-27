utils::globalVariables(c("Year"))

#' Reef Model Parameter Plotting Functions
#'
#' This file contains plotting functions to help users verify and visualize
#' the input parameters of a mizerReef model. These plots are intended for
#' model setup, diagnostics, and parameter checking, not for visualizing
#' simulation results.
#'
#' Functions in this file:
#'   - plotVulnerable: Plots the vulnerability to predation for each species
#'     by weight. Useful for checking refuge and predation profiles.
#'   - plotRefuge: Plots the refuge profile for each species by length,
#'     showing the complement of vulnerability. Useful for verifying refuge
#'     parameterization.
#'   - plotDegScale: Plots a heatmap of degradation scaling parameters for
#'     refuge density across bleaching years and size bins. Useful for
#'     checking input degradation trajectories.
#'
#' Each plot also has an interactive Plotly version for exploring results.
#'
#' These tools are especially useful for:
#'   - Ensuring input data is correctly formatted and interpreted
#'   - Diagnosing unexpected model behavior due to parameter choices
#'   - Comparing parameter sets before running simulations
#'
#' For plots of simulation results and time series, see the summary plotting
#' functions in `reef-summary_plots.R`.
#'
#' @section Usage:
#' Call these functions with a `MizerReefParams` object or relevant parameter
#' data. See individual function documentation for details and examples.
#'
#' @author cmbeese
#' @concept parameterPlots
#' @family plotting functions
#' @seealso [reef-summary_plots.R], [MizerReefParams], [setRefuge()],
#'   [setDegradation()]
NULL

#' Plot the vulnerability to predation of species by weight
#'
#' When called with a \linkS4class{MizerParams} object the initial
#' vulnerability is plotted. The complement of refuge.
#'
#' @param object An object of class \linkS4class{MizerParams}
#'
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a
#'                  numeric vector with the species indices, or a logical
#'                  vector indicating for each species whether it is to be
#'                  selected (TRUE) or not.
#'
#' @param all.sizes If TRUE, then feeding level is plotted also for sizes
#'                  outside a species' size range. Default FALSE.
#'
#' @param return_data   A boolean value that determines whether the formatted
#'                      data used for the plot is returned instead of the plot
#'                      itself. Default value is FALSE.
#' @param ... unused
#'
#' @return  A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'          frame the vulnerability at each size
#'
#' @import ggplot2
#' @import plotly
#' @importFrom dplyr mutate filter select arrange rename summarise
#' @export
#'
#' @family plotting functions
#' @concept refugePlots
#' @seealso [plotting_functions], [setRefuge()], [plotRefuge()]
plotVulnerable <- function(object,
                           species = NULL,
                           all.sizes = FALSE,
                           return_data = FALSE, ...) {
    # input check ----
    assert_that(
        is.flag(all.sizes),
        is.flag(return_data)
    )

    if (is(object, "MizerSim")) {
        ## sim values ----
        params <- object@params
        warning("This functionality is not set up yet.")
    } else if (is(object, "MizerParams")) {
        ## params values ----
        params <- object
    }

    # make plot_dat ----
    sp <- params@species_params

    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }

    # Calculate proportion of fish in refuge
    vul <- getVulnerable(params)

    ## species selector ----
    sel_sp <- valid_species_arg(params, species,
        return.logical = TRUE,
        error_on_empty = TRUE
    )
    species <- dimnames(params@initial_n)$sp[sel_sp]
    species <- gsub("inverts", NA, species)
    species <- species[!is.na(species)]
    sel_sp <- which(!is.na(species))
    vul <- vul[sel_sp, , drop = FALSE]

    ## data frame from selected species -----
    plot_dat <- data.frame(
        w = rep(params@w, each = length(species)),
        value = c(vul),
        Species = species
    )

    if (!all.sizes) {
        # Remove vulnerability for sizes outside a species' size range
        for (sp in species) {
            plot_dat$value[plot_dat$Species == sp &
                (plot_dat$w < params@species_params[sp, "w_min"] |
                    plot_dat$w > params@species_params[sp, "w_max"])] <- NA
        }

        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }

    ## colors ----
    legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
    plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)
    linesize <- rep(0.8, length(legend_levels))
    names(linesize) <- names(params@linetype[legend_levels])

    ## return data if requested ----
    if (return_data) {
        return(plot_dat)
    }

    # plot ----
    p <- ggplot(plot_dat, aes(group = Species))

    ## labels and scales ----
    p + geom_line(aes(
        x = w, y = value,
        colour = Legend, linetype = Legend,
        linewidth = Legend
    )) +
        labs(
            colour = "Species", linetype = "Species",
            linewidth = "Species"
        ) +
        scale_x_continuous(name = "Size [g]", trans = "log10") +
        scale_y_continuous(
            name = "Proportion Vulnerable",
            limits = c(0, 1)
        ) +
        scale_colour_manual(
            values = params@linecolour[legend_levels],
            labels = group_names
        ) +
        scale_linetype_manual(
            values = params@linetype[legend_levels],
            labels = group_names
        ) +
        scale_discrete_manual("linewidth",
            values = linesize,
            labels = group_names
        )
}


#' @rdname plotVulnerable
#' @export
plotlyVulnerable <- function(object,
                             species = NULL,
                             all.sizes = FALSE,
                             return_data = FALSE, ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotVulnerable", argg),
        tooltip = c("Species Group", "w", "value")
    )
}


#' Plot the refuge profile, species by length
#'
#' When called with a \linkS4class{MizerParams} object the initial
#' refuge profile is plotted. The complement of vulnerability.
#'
#' @param object An object of class \linkS4class{MizerParams}
#'
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a
#'                  numeric vector with the species indices, or a logical
#'                  vector indicating for each species whether it is to be
#'                  selected (TRUE) or not.
#'
#' @param all.sizes If TRUE, then feeding level is plotted also for sizes
#'                  outside a species' size range. Default FALSE.
#'
#' @param return_data   A boolean value that determines whether the formatted
#'                      data used for the plot is returned instead of the plot
#'                      itself. Default value is FALSE.
#' @param ... unused
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept refugePlots
#' @seealso [plotting_functions], [setRefuge()], [plotVulnerable()]
plotRefuge <- function(object,
                       species = NULL,
                       all.sizes = FALSE,
                       return_data = FALSE, ...) {
    # input check ----
    assert_that(
        is.flag(all.sizes),
        is.flag(return_data)
    )

    if (is(object, "MizerSim")) {
        ## sim values ----
        params <- object@params
        warning("This functionality is not set up yet.")
    } else if (is(object, "MizerParams")) {
        ## params ----
        params <- object
    }

    # make plot_dat ----
    ## params values ----
    params <- object
    sp <- params@species_params

    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }

    ## Calculate proportion of fish in refuge
    vul <- getVulnerable(params)
    refuge <- (1 - vul)

    ## species selector ----
    sel_sp <- valid_species_arg(params, species,
        return.logical = TRUE,
        error_on_empty = TRUE
    )
    species <- dimnames(params@initial_n)$sp[sel_sp]
    species <- gsub("inverts", NA, species)
    species <- species[!is.na(species)]
    sel_sp <- which(!is.na(species))
    refuge <- refuge[sel_sp, , drop = FALSE]

    # Convert length bins in to weight bins for each functional group
    group_length_bins <- matrix(0,
        nrow = length(species),
        ncol = length(params@w)
    )
    for (i in seq_along(species)) {
        group_length_bins[i, ] <- (params@w / sp$a[i])^(1 / sp$b[i])
    }
    group_length_bins <- group_length_bins[sel_sp, , drop = FALSE]

    # Set x axis limit for plots
    x_limit <- max(sp$l_max)

    ## data frame from selected species -----
    plot_dat <- data.frame(
        w = rep(params@w, each = length(species)),
        l = c(group_length_bins),
        value = c(refuge),
        Species = species
    )

    ## remove sizes outsides range ----
    if (!all.sizes) {
        # Remove vulnerability for sizes outside a species' size range
        for (sp in species) {
            plot_dat$value[plot_dat$Species == sp &
                (plot_dat$w < params@species_params[sp, "w_min"] |
                    plot_dat$w > params@species_params[sp, "w_max"])] <- NA
        }

        plot_dat <- plot_dat[complete.cases(plot_dat), ]
    }

    ## colors ----
    legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
    plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)
    linesize <- rep(0.8, length(legend_levels))
    names(linesize) <- names(params@linetype[legend_levels])

    ## return data if requested ----
    if (return_data) {
        return(plot_dat)
    }

    # plot ----
    ## faceting ----
    p <- ggplot(plot_dat, aes(group = Species)) +
        facet_wrap(~Species, scales = "free_x") +
        theme(
            strip.text.x = element_text(size = 6),
            strip.background = element_rect(fill = "slategray4")
        )

    ## labels and scales ----
    p + geom_line(aes(
        x = l, y = value,
        colour = Legend, linetype = Legend,
        linewidth = Legend
    )) +
        labs(
            colour = "Species", linetype = "Species",
            linewidth = "Species"
        ) +
        scale_x_continuous(
            name = "Total Length (cm)",
            limits = c(0, x_limit)
        ) +
        # scale_x_continuous(name = "Log Size [g]", trans = "log10",
        #                    breaks = c(10^-2, 10^0, 10^2, 10^4),
        #                    labels = c(-2, 0, 2, 4)) +
        scale_y_continuous(name = "Proportion Protected", limits = c(0, 1)) +
        scale_colour_manual(
            values = params@linecolour[legend_levels],
            labels = group_names
        ) +
        scale_linetype_manual(
            values = params@linetype[legend_levels],
            labels = group_names
        ) +
        scale_discrete_manual("linewidth",
            values = linesize,
            labels = group_names
        )
}

#' @rdname plotRefuge
#' @export
plotlyRefuge <- function(object,
                         species = NULL,
                         ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotRefuge", argg),
        tooltip = c("Functional Group", "w", "value")
    )
}

#' Plot heatmap of degradation scaling parameters
#'
#' Creates a heatmap showing the scaling of refuge density across bleaching
#' years and refuge size bins for a given degradation trajectory. The input
#' data should have the first column as bleaching year, and the remaining
#' columns as refuge size bins.
#'
#' @param trajectory Character string: one of "rubble", "algae", "recovery",
#'                   "constant". Determines which built-in trajectory to plot.
#'
#' @param file Optional. Path to CSV file to use. If NULL, uses built-in
#'             data from inst/.
#'
#' @param return_data Logical. If TRUE, returns the formatted data frame
#'                    instead of the plot. Default FALSE.
#'
#' @return ggplot2 object (heatmap), or data frame if return_data = TRUE.
#' @importFrom reshape2 melt
#' @export
#' @family plotting functions
#' @concept parameterPlots
#' @seealso [plotting_functions], [setDegradation()], [reefDegrade()]
plotDegScale <- function(trajectory = NULL,
                         file = NULL,
                         return_data = FALSE) {
    # Validate trajectory
    valid_trajectories <- c("rubble", "algae", "recovery")
    # If trajectory is missing, default to 'rubble'
    if (is.null(trajectory)) {
        trajectory <- "rubble"
    }
    # If trajectory is a vector, use the first value
    if (length(trajectory) > 1) {
        trajectory <- trajectory[1]
    }
    if (!trajectory %in% valid_trajectories) {
        stop(paste0("Invalid trajectory: '", trajectory, "'. Valid options are: ", paste(valid_trajectories, collapse = ", ")))
    }
    # Use built-in R data object if file is NULL, otherwise use CSV
    if (is.null(file)) {
        data_name <- switch(trajectory,
            rubble = "rubble_scale",
            algae = "algae_scale",
            recovery = "recovery_scale"
        )
        # Try to load the data object from the package
        # First, try to get from global env, then from package namespace
        if (!exists(data_name, envir = .GlobalEnv)) {
            data(list = data_name, package = "mizerReef", envir = .GlobalEnv)
        }
        if (exists(data_name, envir = .GlobalEnv)) {
            dat <- get(data_name, envir = .GlobalEnv)
        } else if (exists(data_name, envir = asNamespace("mizerReef"))) {
            dat <- get(data_name, envir = asNamespace("mizerReef"))
        } else {
            stop(paste("Could not find data object:", data_name))
        }
        # If it's a data.frame, convert to matrix
        if (is.data.frame(dat)) {
            dat <- as.matrix(dat)
        }
    } else {
        dat <- as.matrix(read.csv(file, header = FALSE))
        # Check that user-provided data is numeric
        if (!is.numeric(dat)) {
            stop("Provided data must be a numeric matrix or CSV with only numeric values.")
        }
    }
    df_long <- reshape2::melt(dat)
    colnames(df_long) <- c("SizeBin", "Year", "Scaling")
    df_long$SizeBin <- as.factor(df_long$SizeBin)
    # Ensure Year is numeric for labeling
    df_long$Year <- as.numeric(df_long$Year)
    # Optionally return data
    if (return_data) {
        return(df_long)
    }
    # Custom x-axis labels: year 0 = 'Bleach', then 1, 2, ...
    breaks <- sort(unique(df_long$Year))
    labels <- ifelse(breaks == 1, "Bleach", as.character(breaks - 1))
    ggplot(df_long, aes(x = factor(Year, levels = breaks), y = SizeBin, fill = Scaling)) +
        geom_tile(colour = "black", linewidth = 0.25) +
        scale_fill_gradient2(
            name = "Refuge\nDensity\nScaling",
            low = "orangered1",
            mid = "gray90",
            high = "springgreen",
            midpoint = 1,
            limits = c(min(df_long$Scaling, na.rm = TRUE), max(df_long$Scaling, na.rm = TRUE))
        ) +
        scale_x_discrete(
            name = "Years Post Bleaching",
            breaks = breaks,
            labels = labels
        ) +
        labs(
            y = "Refuge Size Bin",
            title = paste(trajectory, "Trajectory"),
            subtitle = "Degradation scaling heatmap"
        ) +
        theme_minimal() +
        theme(
            axis.title.y = element_text(angle = 0, vjust = 0.5),
            legend.key.height = unit(2, "cm")
        )
}

#' Interactive Plotly version of plotDegScale
#'
#' Returns an interactive heatmap for degradation scaling parameters.
#'
#' @inheritParams plotDegScale
#' @return A plotly object
#' @export
plotlyDegScale <- function(trajectory = NULL, file = NULL) {
    ggplotly(plotDegScale(trajectory = trajectory, file = file))
}
