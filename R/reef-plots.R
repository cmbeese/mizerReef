library(ggplot2)
library(plotly)
library(dplyr)

# Set global variables ----
# Not sure if this should be here??
utils::globalVariables(c("Rate", "Consumer", "Producer",
                         "Species", "l", "value", "Legend", "Model",
                         "name1", "name2", "value.y", "value.x",
                         "rel_diff", "scale_bin", "pr"))

#' Plot the biomass of functional groups and unstructured components through
#' time
#'
#' After running a projection, the biomass of each functional group and each
#' unstructured component can be plotted against time. The biomass is
#' calculated within user defined size limits (min_w, max_w, min_l, max_l, see
#' [get_size_range_array()]).
#'
#' @param sim An object of class \linkS4class{MizerSim}
#' 
#' @param species   The groups to be selected. Optional. By default all target
#'                  groups are selected. A vector of groups names, or a numeric 
#'                  vector with the groups indices, or a logical vector 
#'                  indicating for each group whether it is to be selected 
#'                  (TRUE) or not.
#'                  
#' @param y_ticks The approximate number of ticks desired on the y axis.
#' 
#' @param start_time    The first time to be plotted. Default is the beginning 
#'                      of the time series.
#'   
#' @param end_time  The last time to be plotted. Default is the end of the time
#'                  series.
#'                  
#' @param ylim  A numeric vector of length two providing lower and upper limits
#'              for the y axis. Use NA to refer to the existing minimum or 
#'              maximum. Any values below 1e-20 are always cut off.
#'              
#' @param total A boolean value that determines whether the total biomass from
#'              species is plotted as well. Default is FALSE.
#'              
#' @param background    A boolean value that determines whether background 
#'                      species are included. Ignored if the model does not 
#'                      contain background species. Default is TRUE.
#' @param highlight Name or vector of names of the species to be highlighted.
#' 
#' @param return_data   A boolean value that determines whether the formatted 
#'                      data used for the plot is returned instead of the plot 
#'                      itself. Default value is FALSE
#'                      
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

    # User mizer function to create dataframe ----
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

    # Algae ----
    ba <- unlist(sim@n_other$algae)
    dim(ba) <- dim(sim@n_other$algae)
    dimnames(ba) <- dimnames(sim@n_other$algae)
    times <- as.numeric(dimnames(ba)[[1]])
    ba <- ba[(times >= start_time) & (times <= end_time), , drop = FALSE]
    ba <- melt(ba)
    
        # Implement ylim and a minimal cutoff and bring columns in desired order
        min_value <- 1e-20
        ba <- ba[ba$value >= min_value &
                     (is.na(ylim[1]) | ba$value >= ylim[1]) &
                     (is.na(ylim[2]) | ba$value <= ylim[2]), c(1, 3, 2)]
        names(ba) <- c("Year", "Biomass", "Species")
        ba$Legend <- ba$Species
    
    # Detritus ----
    bd <- unlist(sim@n_other$detritus)
    dim(bd) <- dim(sim@n_other$detritus)
    dimnames(bd) <- dimnames(sim@n_other$detritus)
    times <- as.numeric(dimnames(bd)[[1]])
    bd <- bd[(times >= start_time) & (times <= end_time), , drop = FALSE]
    bd <- melt(bd)

        # Implement ylim and a minimal cutoff and bring columns in desired order
        min_value <- 1e-20
        bd <- bd[bd$value >= min_value &
                     (is.na(ylim[1]) | bd$value >= ylim[1]) &
                     (is.na(ylim[2]) | bd$value <= ylim[2]), c(1, 3, 2)]
        names(bd) <- c("Year", "Biomass", "Species")
        bd$Legend <- bd$Species

    plot_dat <- rbind(df, ba, bd)
    if (return_data) return(plot_dat)

    plotDataFrame(plot_dat, params, xlab = "Year", ylab = "Biomass [g]",
                  ytrans = "log10",
                  y_ticks = y_ticks, highlight = highlight,
                  legend_var = "Functional Group")
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
             tooltip = c("Functional Group", "Year", "Biomass"))
}

#' Plot the vulnerability to predation of species by size
#'
#' When called with a \linkS4class{MizerParams}
#' object the initial vulnerability is plotted.
#'
#' @param object An object of class \linkS4class{MizerParams}
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a 
#'                  numeric vector with the species indices, or a logical 
#'                  vector indicating for each species whether it is to be 
#'                  selected (TRUE) or not.
#' @param all.sizes If TRUE, then feeding level is plotted also for sizes
#'                  outside a species' size range. Default FALSE.
#' @param return_data   A boolean value that determines whether the formatted 
#'                      data used for the plot is returned instead of the plot 
#'                      itself. Default value is FALSE.
#' @param ... unused
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept plots
#' @seealso [plotting_functions], [setRefuge()]
plotVulnerable <- function(object, 
                           species = NULL,
                           all.sizes = FALSE,
                           return_data = FALSE,...) {
    
    # input check ----
    assert_that(is.flag(all.sizes),
                is.flag(return_data))
    
    if (is(object, "MizerSim")) {
        ## sim values ----
        params <- object@params
        warning('This functionality is not set up yet.')
    } else if (is(object, "MizerParams")) {
        ## params values ----
        params <- object
    }
    
    # make plot_dat ----
    sp <- params@species_params
    no_sp <- dim(params@interaction)[1]
    
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }
    
    # Calculate proportion of fish in refuge
    vul <- getVulnerable(params)
    refuge <- (1-vul)
    
    ## species selector ----
    sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                error_on_empty = TRUE)
    species <- dimnames(params@initial_n)$sp[sel_sp]
    species <- gsub('inverts',NA,species)
    species <- species[!is.na(species)]
    sel_sp <- which(!is.na(species))
    refuge <- refuge[sel_sp, , drop = FALSE]
    
    # Convert length bins in to weight bins for each functional group
    group_length_bins <- matrix(0, nrow = length(species), ncol = length(params@w))
    for (i in 1:length(species)) {
        group_length_bins[i,] <- (params@w / sp$a[i])^(1 / sp$b[i])
    }
    group_length_bins <- group_length_bins[sel_sp, , drop = FALSE]
    
    # Set x axis limit for plots
    x_limit = max(sp$l_max)
    
    ## data frame from selected species -----
    plot_dat <- data.frame(w = rep(params@w, each = length(species)),
                           l = c(group_length_bins),
                           value = c(refuge),
                           Species = species)
    
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
    if (return_data) return(plot_dat)
    
    # plot ----
    ## faceting ----
    p <- ggplot(plot_dat, aes(group = Species)) +
        facet_wrap(~ Species, scales = "free_x") +
        theme(strip.text.x = element_text(size = 6))
    #   strip.background = element_blank(),
    #   strip.text.x =element_blank())
    
    ## labels and scales ----
    p + geom_line(aes(x = l, y = value,
                      colour = Legend, linetype = Legend,
                      linewidth = Legend)) +
        labs(colour = 'Functional Group', linetype = 'Functional Group',
             linewidth = 'Functional Group') +
        scale_x_continuous(name = "Total Length [cm]",
                           limits = c(0,x_limit)) +
        # scale_x_continuous(name = "Log Size [g]", trans = "log10",
        #                    breaks = c(10^-2, 10^0, 10^2, 10^4),
        #                    labels = c(-2, 0, 2, 4)) +
        scale_y_continuous(name = "Proportion Protected", limits = c(0, 1)) +
        scale_colour_manual(values = params@linecolour[legend_levels],
                            labels = group_names) +
        scale_linetype_manual(values = params@linetype[legend_levels],
                              labels = group_names) +
        scale_discrete_manual("linewidth", values = linesize,
                              labels = group_names)
}


#' @rdname plotVulnerable
#' @export
plotlyVulnerable <- function(object,
                             species = NULL,
                             all.sizes = FALSE,
                             return_data = FALSE,...) {

    argg <- as.list(environment())
    ggplotly(do.call("plotVulnerable", argg),
             tooltip = c("Functional Group", "w", "value"))
}


#' Plot the refuge profile, species by size
#'
#' When called with a \linkS4class{MizerParams}
#' object the initial vulnerability is plotted.
#'
#' @param object An object of class \linkS4class{MizerParams}
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a 
#'                  numeric vector with the species indices, or a logical 
#'                  vector indicating for each species whether it is to be 
#'                  selected (TRUE) or not.
#' @param all.sizes If TRUE, then feeding level is plotted also for sizes
#'                  outside a species' size range. Default FALSE.
#' @param return_data   A boolean value that determines whether the formatted 
#'                      data used for the plot is returned instead of the plot 
#'                      itself. Default value is FALSE.
#' @param ... unused
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept plots
#' @seealso [plotting_functions], [setRefuge()]
plotRefuge <- function(object, 
                       species = NULL,
                       all.sizes = FALSE,
                       return_data = FALSE,...) {
    
    # input check ----
    assert_that(is.flag(all.sizes),
                is.flag(return_data))
    
    if (is(object, "MizerSim")) {
        ## sim values ----
        params <- object@params
        warning('This functionality is not set up yet.')
    } else if (is(object, "MizerParams")) {
        ## params ----
        params <- object
    }
    
    # make plot_dat ----
        ## params values ----
        params <- object
        sp <- params@species_params
        no_sp <- dim(params@interaction)[1]
        
        if (is.null(params@species_params$group_names)){
            group_names <- params@species_params$species
        } else {
            group_names <- params@species_params$group_names
        }
        
        ## Calculate proportion of fish in refuge
        vul <- getVulnerable(params)
        refuge <- (1-vul)

        ## species selector ----
        sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                    error_on_empty = TRUE)
        species <- dimnames(params@initial_n)$sp[sel_sp]
        species <- gsub('inverts',NA,species)
        species <- species[!is.na(species)]
        sel_sp <- which(!is.na(species))
        refuge <- refuge[sel_sp, , drop = FALSE]
        
        # Convert length bins in to weight bins for each functional group
        group_length_bins <- matrix(0, nrow = length(species), ncol = length(params@w))
        for (i in 1:length(species)) {
            group_length_bins[i,] <- (params@w / sp$a[i])^(1 / sp$b[i])
        }
        group_length_bins <- group_length_bins[sel_sp, , drop = FALSE]
        
        # Set x axis limit for plots
        x_limit = max(sp$l_max)
        
        ## data frame from selected species -----
        plot_dat <- data.frame(w = rep(params@w, each = length(species)),
                               l = c(group_length_bins),
                               value = c(refuge),
                               Species = species)
        
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
        if (return_data) return(plot_dat)
        
        # plot ----
        ## faceting ----
        p <- ggplot(plot_dat, aes(group = Species)) +
            facet_wrap(~ Species, scales = "free_x") +
            theme(strip.text.x = element_text(size = 6))
            #   strip.background = element_blank(),
            #   strip.text.x =element_blank())
        
        ## labels and scales ----
        p + geom_line(aes(x = l, y = value,
                          colour = Legend, linetype = Legend,
                          linewidth = Legend)) +
            labs(colour = 'Functional Group', linetype = 'Functional Group',
                 linewidth = 'Functional Group') +
            scale_x_continuous(name = "Total Length (cm)",
                               limits = c(0,x_limit)) +
            # scale_x_continuous(name = "Log Size [g]", trans = "log10",
            #                    breaks = c(10^-2, 10^0, 10^2, 10^4),
            #                    labels = c(-2, 0, 2, 4)) +
            scale_y_continuous(name = "Proportion Protected", limits = c(0, 1)) +
            scale_colour_manual(values = params@linecolour[legend_levels],
                                labels = group_names) +
            scale_linetype_manual(values = params@linetype[legend_levels],
                                  labels = group_names) +
            scale_discrete_manual("linewidth", values = linesize,
                                  labels = group_names)
}

#' @rdname plotRefuge
#' @export
plotlyRefuge <- function(object,
                             species = NULL,
                             ...) {

    argg <- as.list(environment())
    ggplotly(do.call("plotRefuge", argg),
             tooltip = c("Functional Group", "w", "value"))
}


#' Plot the total productivity for each functional group
#'
#' When called with a \linkS4class{MizerParams} object the steady state 
#' productivity is plotted. When called with a \linkS4class{MizerSim} 
#' object the productivity is plotted.
#'
#' @param object An object of class \linkS4class{MizerParams}
#' 
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a 
#'                  numeric vector with the species indices, or a logical 
#'                  vector indicating for each species whether it is to be 
#'                  selected (TRUE) or not.
#'                  
#' @param start_time    The first time to be plotted. Default is the beginning 
#'                      of the time series.
#'                      
#' @param end_time  The last time to be plotted. Default is the end of the time
#'                  series.
#'                  
#' @param return_data   A boolean value that determines whether the formatted 
#'                      data used for the plot is returned instead of the plot 
#'                      itself. Default value is FALSE.
#'                      
#' @param min_fishing_l parameters be passed to [getProductivity()]. The 
#'                      minimum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to 7 cm.
#'                      
#' @param max_fishing_l parameters be passed to [getProductivity()]. The 
#'                      maximum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to max length. 
#' @param ... unused
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept plots
#' @seealso [plotting_functions], [setRefuge()]
plotProductivity <- function(object,
                             species = NULL,
                             return_data = FALSE,
                             min_fishing_l = NULL, max_fishing_l = NULL,
                             start_time = NULL, end_time = NULL,...) {

    if (is(object, "MizerSim")) {
        # sim values ----
        params <- object@params
        warning('This functionality is not set up yet.')
    }
    if (is(object, "MizerParams")) {
        # params ----
        params <- object
    } else {
        stop('object should be a mizerParams or mizerSim object.')
    }
    
    # create plot_dat ----
    
    ## values from object ----
    params <- object
    sp <- params@species_params
    no_sp <- dim(params@interaction)[1]
    
    ## group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }
    
    ## get productivity ----
    prod <- getProductivity(params, 
                            min_fishing_l = min_fishing_l,
                            max_fishing_l = max_fishing_l)
    
    ## species selector ----
    sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                error_on_empty = TRUE)
    species <- dimnames(params@initial_n)$sp[sel_sp]
    species <- gsub('inverts', NA, species)
    species <- species[!is.na(species)]
    sel_sp <- which(!is.na(species))
    prod <- prod[sel_sp, drop = FALSE]
    
    ## data frame from selected species ----
    plot_dat <- data.frame(value = prod, Species = species)

    ## colors ----
    legend_levels   <- intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
    
    ## return data if requested ----
    if (return_data) return(plot_dat)
    
    # plot ----
    p <- ggplot(plot_dat, aes(x = Species, y = value,
                              group = Legend, fill = Legend))
    p + geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(name = expression("Productivity (g/m^2/year)")) +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names)
    
}

#' @rdname plotProductivity
#' @export
plotlyProductivity <- function(object,
                               species = NULL,...) {

    argg <- as.list(environment())
    ggplotly(do.call("plot2Productivity", argg),
             tooltip = c("Species", "value"))
}


#' Plot the productivity of two models in the same plot
#'
#' When called with a \linkS4class{MizerParams}
#' object the steady state productivities are plotted.
#'
#' @param object1 First MizerParams or MizerSim object.
#' 
#' @param object2 Second MizerParams or MizerSim object.
#' 
#' @param name1 An optional string with the name for the first model, to be 
#'              used in the legend. Set to "First" by default.
#'              
#' @param name2 An optional string with the name for the second model, to be
#'              used in the legend. Set to "Second" by default.
#'              
#' @param min_fishing_l parameters be passed to [getProductivity()]. The 
#'                      minimum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to 7 cm.
#'                      
#' @param max_fishing_l parameters be passed to [getProductivity()]. The 
#'                      maximum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to max length.
#'                       
#' @param ... Parameters passed to `plotProductivity()`
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept plots
#' @seealso [plotting_functions], [setRefuge()]
plot2Productivity <- function(object1, object2, 
                              name1 = "First", name2 = "Second",
                              min_fishing_l = NULL, max_fishing_l = NULL,...) {
    
    # get data frames with plotProductivity ----
    sf1 <- plotProductivity(object1, return_data = TRUE, ...)
    sf1$Model <- name1
    sf2 <- plotProductivity(object2, return_data = TRUE, ...)
    sf2$Model <- name2
    sf <- rbind(sf1, sf2)
    
    # if sim, get params ----
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    # group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }
    
    # plot ----
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    
    p <- ggplot(sf, aes(x = Species, y = value, 
                        group = Model, alpha = Model,
                        fill = Legend))
    p + geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_y_continuous(name = expression("Productivity (g/m^2/year)")) +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names) +
        labs(fill = "Functional Groups") + 
        scale_alpha_manual(values = c(0.5,1),
                           labels = c(name1, name2))
    
}

#' @rdname plot2Productivity
#' @export
plotly2Productivity <- function(object,
                                species = NULL,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("plot2Productivity", argg),
             tooltip = c("Species", "value"))
}

#' Plot the relative difference between the productivity rates of two models 
#' 
#' This function creates a barplot with the difference between the productivity 
#' rates for each functional group in two models. Usually used in conjunction
#' with [plotTotalBiomassRelative()] to check for decoupling.
#' 
#' The individual productivity rates are calculated by the 
#' [plotProductivity()] function which is passed all additional arguments you 
#' supply. So you can for example set a size range with the `min_fishing_l`
#' and `max_fishing_l` arguments. See [plotProductivity()] for more details.
#' 
#' @param object1 An object of class MizerSim or MizerParams
#' 
#' @param object2 An object of class MizerSim or MizerParams
#' 
#' @param return_data   A boolean value that determines whether the formatted 
#'                      data used for the plot is returned instead of the plot 
#'                      itself. Default value is FALSE.
#' 
#' @param ... Parameters passed to `plotSpectra()`
#' 
#' @return A ggplot2 object
#' @concept plots
#' @export
plotProductivityRelative <- function(object1, object2,
                                     return_data = FALSE,...) {
    
    # get data frames with plotProductivity ----
    sf1 <- plotProductivity(object1, return_data = TRUE, ...)
    sf2 <- plotProductivity(object2, return_data = TRUE, ...)
    
    # Calclulate relative difference
    sf <- dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |>
          dplyr::mutate(rel_diff = (value.y - value.x) / (value.x + value.y))
    
    if(return_data == TRUE){return(sf)}
    
    # object check ----
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    # group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }
    
    # plot -----
    legend_levels <- intersect(names(params@linecolour),
                               unique(sf$Legend))
    
    p <- ggplot(sf, aes(x = Species, y = rel_diff,
                        fill = Legend))
    p + geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_y_continuous(name = "Relative Difference") +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names) +
        labs(fill = "Functional Groups") + 
        geom_hline(yintercept = 0, linetype = 1,
                   colour = "dark grey", linewidth = 0.75)
}

#' @rdname plotProductivityRelative
#' @export
plotlyProductivityRelative <- function(object,
                                species = NULL,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("plotProductivityRelative", argg),
             tooltip = c("Species", "value"))
}

#' Plot the total fishable biomass for each functional group at steady state
#' 
#' This functions creates a barplot with the biomass of each functional group
#' within a size range. Usually used in conjunction with [plotProductivity()] 
#' to check for decoupling.  
#'
#' @param object An object of class \linkS4class{MizerParams}
#' 
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a 
#'                  numeric vector with the species indices, or a logical 
#'                  vector indicating for each species whether it is to be 
#'                  selected (TRUE) or not.
#'                  
#' @param min_fishing_l parameters be passed to [getProductivity()]. The 
#'                      minimum length (cm) of fished individuals for
#'                      biomass estimates. Defaults to 7 cm.
#'                      
#' @param max_fishing_l parameters be passed to [getProductivity()]. The 
#'                      maximum length (cm) of fished individuals for
#'                      biomass estimates. Defaults to max length.
#'                  
#' @param return_data   A boolean value that determines whether the formatted 
#'                      data used for the plot is returned instead of the plot 
#'                      itself. Default value is FALSE.
#'                      
#' @param ... unused
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept plots
#' @seealso [plotting_functions], [plotBiomass()]
plotTotalBiomass <- function(object,
                             species = NULL,
                             min_fishing_l = NULL, max_fishing_l = NULL,
                             return_data = FALSE, ...) {
    
    # object checks ----
    if (is(object, "MizerSim")) {
        ## sim values ----
        params <- object@params
        warning('This functionality is not set up yet.')
    }
    if (is(object, "MizerParams")) {
        ## params ----
        params <- object
    } else {
        stop('object should be a mizerParams or mizerSim object.')
    }
    
    # create plot_dat ----
    
    ## values from object ----
    params <- object
    sp <- params@species_params
    no_sp <- dim(params@interaction)[1]
    
    ## group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }
    
    ## get productivity ----
    biom <- getBiomass(params,
                       min_l = min_fishing_l,
                       max_l = max_fishing_l)
    
    ## species selector ----
    sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                error_on_empty = TRUE)
    species <- dimnames(params@initial_n)$sp[sel_sp]
    species <- gsub('inverts', NA, species)
    species <- species[!is.na(species)]
    sel_sp <- which(!is.na(species))
    biom <- biom[sel_sp, drop = FALSE]
    
    ## data frame from selected species ----
    plot_dat <- data.frame(value = biom, Species = species)
    
    ## colors ----
    legend_levels   <- intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
    
    ## return data if requested ----
    if (return_data) return(plot_dat)
    
    # plot ----
    p <- ggplot(plot_dat, aes(x = Species, y = value,
                              group = Legend, fill = Legend))
    p + geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(name = expression("Productivity (g/m^2/year)")) +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names)

}

#' @rdname plotTotalBiomass
#' @export
plotlyTotalBiomass <- function(object,
                               species = NULL,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("plotTotalBiomass", argg),
             tooltip = c("Species", "value"))
}


#' Plot the total biomass of two models in the same plot
#'
#' When called with a \linkS4class{MizerParams}
#' object the steady state biomasses are plotted.
#' 
#' @param object1 First MizerParams or MizerSim object.
#' 
#' @param object2 Second MizerParams or MizerSim object.
#' 
#' @param name1 An optional string with the name for the first model, to be 
#'              used in the legend. Set to "First" by default.
#'              
#' @param name2 An optional string with the name for the second model, to be
#'              used in the legend. Set to "Second" by default.
#'              
#' @param min_fishing_l parameters be passed to [getProductivity()]. The 
#'                      minimum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to 7 cm.
#'                      
#' @param max_fishing_l parameters be passed to [getProductivity()]. The 
#'                      maximum length (cm) of fished individuals for
#'                      productivity estimates. Defaults to max length.
#'                       
#' @param ... Parameters to pass to [plotProductivity()] or [plotTotalBiomass()]
#'
#' @return A ggplot2 object
#'
#' @export
#' @family plotting functions
#' @concept plots
#' @seealso [plotting_functions], [setRefuge()]
plot2TotalBiomass <- function(object1, object2, 
                              name1 = "First", name2 = "Second",
                              min_fishing_l = NULL, max_fishing_l = NULL,...) {
    
    # get data frames with plotProductivity ----
    sf1 <- plotTotalBiomass(object1, return_data = TRUE, ...)
    sf1$Model <- name1
    sf2 <- plotTotalBiomass(object2, return_data = TRUE, ...)
    sf2$Model <- name2
    sf <- rbind(sf1, sf2)
    
    # if sim, get params ----
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    # group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }
    
    # plot ----
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    
    p <- ggplot(sf, aes(x = Species, y = value, 
                        group = Model, alpha = Model,
                        fill = Legend))
    p + geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_y_continuous(name = expression("Biomass (g/m^2)")) +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names) +
        labs(fill = "Functional Groups") + 
        scale_alpha_manual(values = c(0.5,1),
                           labels = c(name1, name2))
    
}

#' @rdname plot2TotalBiomass
#' @export
plotly2TotalBiomass <- function(object,
                               species = NULL,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("plot2TotalBiomass", argg),
             tooltip = c("Species", "value"))
}

#' Plot the relative difference in between the total fishable biomasses of each
#' each functional group at steady state
#' 
#' This functions creates a barplot with the relative change in biomass of 
#' each functional group within a size range. Usually used in conjunction 
#' with [plotProductivityRelative()] to check for decoupling. 
#' 
#' The individual biomasses are calculated by the [plotTotalBiomass()] function 
#' which is passed all additional arguments you supply. So you can for example 
#' set a size range with the `min_fishing_l` and `max_fishing_l` arguments. 
#' See [plotTotalBiomass()] for more details.
#' 
#' @param object1 An object of class MizerSim or MizerParams
#' @param object2 An object of class MizerSim or MizerParams
#' @param ... Parameters passed to `plotSpectra()`
#' 
#' @return A ggplot2 object
#' @concept plots
#' @export
plotTotalBiomassRelative <- function(object1, object2, ...) {

    # get data frames with plotProductivity ----
    sf1 <- plotTotalBiomass(object1, return_data = TRUE, ...)
    sf2 <- plotTotalBiomass(object2, return_data = TRUE, ...)
    
    # Calculate relative difference
    sf <- dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |>
          dplyr::mutate(rel_diff = (value.y - value.x) / (value.x + value.y))
    
    # if sim, get params ----
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    # group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
    }
    
    # plot ----
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    
    p <- ggplot(sf, aes(x = Species, y = rel_diff, fill = Legend))
    
    p + geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_y_continuous(name = "Relative Difference") +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names) +
        labs(fill = "Functional Groups") + 
        geom_hline(yintercept = 0, linetype = 1,
                   colour = "dark grey", linewidth = 0.75)
    
}
    
#' @rdname plotTotalBiomassRelative
#' @export
plotlyTotalBiomassRelative <- function(object,
                                species = NULL,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("TotalBiomassRelative", argg),
             tooltip = c("Species", "value"))
}