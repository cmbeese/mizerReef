library(ggplot2)
library(plotly)
library(dplyr)

#' Plot the vulnerability to predation of species by size
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
#' @export
#' 
#' @family plotting functions
#' @concept refugePlots
#' @seealso [plotting_functions], [setRefuge()], [plotRefuge()]
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
    
    # Set x axis limit for plots
    x_limit = max(sp$l_max)
    
    ## data frame from selected species -----
    plot_dat <- data.frame(w = rep(params@w, each = length(species)),
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
    p <- ggplot(plot_dat, aes(group = Species))
    #    facet_wrap(~ Species, scales = "free_x") +
    #    theme(strip.text.x = element_text(size = 6))
    #   strip.background = element_blank(),
    #   strip.text.x =element_blank())
    
    ## labels and scales ----
    p + geom_line(aes(x = w, y = value,
                      colour = Legend, linetype = Legend,
                      linewidth = Legend)) +
        labs(colour = 'Functional Group', linetype = 'Functional Group',
             linewidth = 'Functional Group') +
        # scale_x_continuous(name = "Total Length [cm]",
        #                    limits = c(0,x_limit)) +
        scale_x_continuous(name = "Log Size [g]", trans = "log10") +#,
                           #breaks = c(10^-2, 10^0, 10^2, 10^4),
                           #labels = c(-2, 0, 2, 4)) +
        scale_y_continuous(name = "Proportion Protected", 
                           limits = c(0, 1)) +
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