utils::globalVariables(c("Year", "SizeBin"))

#' Plot the vulnerability to predation of species by weight
#'
#' When called with a \linkS4class{MizerParams} object the initial
#' vulnerability is plotted. The complement of refuge.
#'
#' @param object An object of class \linkS4class{MizerParams} or
#'   \linkS4class{MizerSim}. If a \linkS4class{MizerSim} object is provided, the
#'   abundance at the last time step is used.
#'
#' @param species The species to be selected. Optional. By default all species
#'   are selected. A vector of species names, or a numeric vector with the
#'   species indices, or a logical vector indicating for each species whether
#'   it is to be selected (TRUE) or not.
#'
#' @param all.sizes If TRUE, then vulnerability is plotted also for sizes
#'   outside a species' size range. Default FALSE.
#'
#' @param time_step If `object` is a \linkS4class{MizerSim} object, this
#'   optional parameter specifies which time step to use for calculating
#'   vulnerability. Default is the last time step.
#'
#' @param return_data A boolean value that determines whether the formatted
#'   data used for the plot is returned instead of the plot itself. Default
#'   value is FALSE.
#'
#' @param ... unused
#'
#' @return A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'   frame the vulnerability at each size
#'
#' @import ggplot2
#' @importFrom dplyr mutate filter select arrange rename summarise
#' @export
#'
#' @title Plot the vulnerability to predation of species by weight
#' @name plotVulnerable
#' @family plotting functions
#' @concept refugePlots
#' @seealso [plotting_functions], [setRefuge()], [plotRefugeProfile()]
plotVulnerable <- function(object,
                           species = NULL,
                           all.sizes = FALSE,
                           time_step = NULL,
                           return_data = FALSE, ...) {
    # input check ----
    assert_that(
        is.flag(all.sizes),
        is.flag(return_data)
    )

    if (is(object, "MizerSim")) {
        ## sim values ----
        params <- object@params
        # Get time for end of simulation
        if (is.null(time_step)) {
            t <- max(as.numeric(dimnames(object@n)$time))
        }
        t <- time_step
        vul <- getVulnerable(object, time_range = t)
    } else if (is(object, "MizerParams")) {
        ## params values ----
        params <- object

        # Calculate proportion of fish in refuge
        vul <- getVulnerable(params)
    }

    # make plot_dat ----
    sp <- params@species_params

    if (is.null(params@species_params$group_names)) {
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }

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


#' @importFrom plotly ggplotly
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
#' When called with a \linkS4class{MizerParams} object the initial refuge
#' profile is plotted. The complement of vulnerability.
#'
#' @param object An object of class \linkS4class{MizerParams}
#'
#' @param species The species to be selected. Optional. By default all species
#'   are selected. A vector of species names, or a numeric vector with the
#'   species indices, or a logical vector indicating for each species whether
#'   it is to be selected (TRUE) or not.
#'
#' @param all.sizes If TRUE, then feeding level is plotted also for sizes
#'   outside a species' size range. Default FALSE.
#'
#' @param return_data A boolean value that determines whether the formatted
#'   data used for the plot is returned instead of the plot itself. Default
#'   value is FALSE.
#' @param ... unused
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept refugePlots
#' @seealso [plotting_functions], [setRefuge()], [plotVulnerable()]
plotRefugeProfile <- function(object,
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
        warning("You are plotting the refuge profile from the steady state. To view refuge density through time, use plotRefugeDensity.")
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

#' @importFrom plotly ggplotly
#' @rdname plotRefugeProfile
#' @export
plotlyRefugeProfile <- function(object,
                                species = NULL,
                                ...) {
    argg <- as.list(environment())
    ggplotly(do.call("plotRefugeProfile", argg),
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
#' This function is flexible in its input:
#' \itemize{
#'   \item If \code{trajectory} is provided, it takes precedence. It can be:
#'     \itemize{
#'       \item A character string (\code{"rubble"}, \code{"algae"},
#'         \code{"recovery"}) to use built-in data.
#'       \item A user-provided numeric matrix or data.frame (bleaching year in
#'         first column, remaining columns as refuge size bins).
#'     }
#'   \item If \code{trajectory} is not provided, \code{object} must be
#'     supplied and can be a \linkS4class{MizerParams} or
#'     \linkS4class{MizerSim} object. The function will extract the
#'     degradation scaling or trajectory from
#'     \code{object@other_params$refuge_params$trajectory} or
#'     \code{object@other_params$refuge_params$deg_scale}.
#' }
#'
#' @param object An optional object of class \linkS4class{MizerParams} or
#'   \linkS4class{MizerSim}. If provided, the function will attempt to extract
#'   the degradation scaling or trajectory from its \code{@other_params$refuge_params} slot.
#'
#' @param trajectory Optional. Either a character string (\code{"rubble"},
#'   \code{"algae"}, \code{"recovery"}) to use built-in data, or a
#'   user-provided numeric matrix/data.frame with the correct format (bleaching
#'   year in first column, remaining columns as refuge size bins).
#'
#' @param return_data Logical. If TRUE, returns the formatted data frame
#'   instead of the plot. Default FALSE.
#'
#' @return A ggplot2 object (heatmap), or a data frame if
#'   \code{return_data = TRUE}.
#' @importFrom reshape2 melt
#' @export
#' @family plotting functions
#' @concept refugePlots
#' @seealso [plotting_functions], [setDegradation()], [reefDegrade()]
plotDegradationScale <- function(object = NULL,
                                 trajectory = NULL,
                                 return_data = FALSE) {
    # At least a MizerParams object, MizerSim object, or custom trajectory must be provided
    if (is.null(object) && is.null(trajectory)) {
        stop("Either 'object' or 'trajectory' must be provided.")
    }
    valid_trajectories <- c("rubble", "algae", "recovery")
    dat <- NULL

    # 1. If trajectory is provided, it takes precedence
    if (!is.null(trajectory)) {
        # If trajectory is a matrix or data.frame, use it directly
        if (is.matrix(trajectory) || is.data.frame(trajectory)) {
            dat <- as.matrix(trajectory)
            if (!is.numeric(dat)) {
                stop("Provided trajectory matrix/data.frame must be numeric.")
            }
        } else if (is.character(trajectory)) {
            # If trajectory is a vector, use the first value
            if (length(trajectory) > 1) {
                trajectory <- trajectory[1]
            }
            if (!trajectory %in% valid_trajectories) {
                stop(paste0("Invalid trajectory: '", trajectory, "'. Valid options are: ", paste(valid_trajectories, collapse = ", ")))
            }
            data_name <- switch(trajectory,
                rubble = "rubble_scale",
                algae = "algae_scale",
                recovery = "recovery_scale"
            )
            # Try to load the data object from the package
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
            if (is.data.frame(dat)) {
                dat <- as.matrix(dat)
            }
        } else {
            stop("'trajectory' must be a character, matrix, or data.frame.")
        }
    } else if (!is.null(object)) {
        # 2. If only object is provided
        if (inherits(object, "MizerSim")) {
            params <- object@params
        } else if (inherits(object, "MizerParams")) {
            params <- object
        } else {
            stop("'object' must be a MizerParams or MizerSim object.")
        }
        # Try to extract trajectory or deg_scale from params@other_params$refuge_params
        if (!is.null(params@other_params$refuge_params$trajectory)) {
            dat <- params@other_params$refuge_params$trajectory
        } else if (!is.null(params@other_params$refuge_params$deg_scale)) {
            dat <- params@other_params$refuge_params$deg_scale
        } else {
            stop("No degradation scaling or trajectory found in object@other_params$refuge_params.")
        }
    }

    # Validate dat
    if (is.null(dat) || !is.matrix(dat) || !is.numeric(dat)) {
        stop("Could not extract a valid numeric matrix for degradation scaling.")
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
        scale_fill_gradientn(
            name = "Refuge\nDensity\nScaling",
            colors = c("deeppink3", "orangered1", "gray90", "springgreen", "dodgerblue2"),
            values = scales::rescale(c(0, 0.95, 1, 1.05, 2)),
            limits = c(0, 2),
            labels = function(x) paste0(round((x - 1) * 100), "%")
        ) +
        scale_x_discrete(
            name = "Years Post Bleaching",
            breaks = breaks,
            labels = labels
        ) +
        scale_y_discrete(
            name = "Refuge Size Bin",
            limits = rev(levels(df_long$SizeBin))
        ) +
        coord_fixed(ratio = 1) +
        labs(
            title = if (is.character(trajectory) && trajectory %in% valid_trajectories) paste(trajectory, "trajectory") else "Custom Trajectory",
        ) +
        theme_minimal() +
        theme(
            axis.title.y = element_text(angle = 0, vjust = 0.5),
            legend.key.height = unit(2, "cm")
        )
}
#' Interactive Plotly version of plotDegradationScale
#'
#' Returns an interactive heatmap for degradation scaling parameters.
#'
#' @inheritParams plotDegradationScale
#' @return A plotly object
#' @importFrom plotly ggplotly
#' @export
plotlyDegradationScale <- function(object = NULL, trajectory = NULL, return_data = FALSE, ...) {
    ggplotly(plotDegradationScale(object = object, trajectory = trajectory, return_data = return_data, ...))
}
