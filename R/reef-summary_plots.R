library(ggplot2)
library(plotly)
#library(dplyr)

# Set global variables ----
# Variables used - ggplot2 - known bug
# Hackiness to get past the 'no visible binding ... ' warning when running check
utils::globalVariables(c("Species", "value", "Model", "Legend",
                         "value.y", "value.x","rel_diff", "l",
                         "y_ticks", "highlight", "Metric",
                         # variables used by ggplot for detritus and algae
                         "Rate", "Source", "Consumer"))

#' Plot the biomass of Species Groups and unstructured components through
#' time
#'
#' After running a projection, the biomass of each Species Group and each
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
#' @return  A ggplot2 object, unless `return_data = TRUE`, in which case a data
#'          frame with the four variables 'Year', 'Biomass', 'Species', 
#'          'Legend' is returned.
#'
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
#' 
#' @concept sumplots
#' @family plotting functions
#' 
#' @seealso [plotBiomass()], [plot2TotalBiomass()], [plotTotalBiomassRelative()],
#'          [plotProductivity()], [plot2Productivity()], [plotProductivityRelative()]
plotBiomass <- function(sim, species = NULL,
                        start_time, end_time,
                        y_ticks = 6, ylim = c(NA, NA),
                        total = FALSE, background = TRUE,
                        highlight = NULL, return_data = FALSE,
                        ...) {
    
    params <- sim@params
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
    
    species <- valid_species_arg(sim, species)
    if (missing(start_time)) start_time <-
        as.numeric(dimnames(sim@n)[[1]][1])
    if (missing(end_time)) end_time <-
        as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]])
    if (start_time >= end_time) {
        stop("start_time must be less than end_time")
    }
    
    # other components
    bu <- unlist(sim@n_other)
    dim(bu) <- dim(sim@n_other)
    dimnames(bu) <- dimnames(sim@n_other)
    times <- as.numeric(dimnames(bu)[[1]])
    bu <- bu[(times >= start_time) & (times <= end_time), , drop = FALSE]
    bu <- reshape2::melt(bu)
    # Implement ylim and a minimal cutoff and bring columns in desired order
    min_value <- 1e-20
    bu <- bu[bu$value >= min_value &
                 (is.na(ylim[1]) | bu$value >= ylim[1]) &
                 (is.na(ylim[2]) | bu$value <= ylim[2]), c(1, 3, 2)]
    names(bu) <- c("Year", "Biomass", "Species")
    bu$Legend <- bu$Species
    
    # Return data ----
    plot_dat <- rbind(df, bu)
    if (return_data) return(plot_dat)
    
    p <- mizer::plotDataFrame(plot_dat, params,
                              xlab = "Year", ylab = "Biomass (g)",
                              ytrans = "log10", y_ticks = y_ticks, 
                              highlight = highlight,
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

#' Show two size spectra in the same plot
#' 
#' @param object1 First MizerParams or MizerSim object.
#' @param object2 Second MizerParams or MizerSim object.
#' @param name1 An optional string with the name for the first model, to be used
#'   in the legend. Set to "First" by default.
#' @param name2 An optional string with the name for the second model, to be
#'   used in the legend. Set to "Second" by default.
#' @param power The abundance is plotted as the number density times the weight
#'   raised to this power. The default power = 1 gives the biomass density,
#'   whereas power = 2 gives the biomass density with respect to logarithmic
#'   size bins.
#' @param ... Parameters to pass to `plotSpectra()`
#' @return A ggplot2 object
#' @export
#' @concept sumplots
plotSpectra2 <- function(object1, object2, 
                         name1 = "First", name2 = "Second",
                         power = 1, ...) {
    
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    sf1 <- plotSpectra(object1, 
                       power = power, return_data = TRUE, ...)
    sf1$Model <- name1
    sf2 <- plotSpectra(object2, 
                       power = power, return_data = TRUE, ...)
    sf2$Model <- name2
    sf <- rbind(sf1, sf2)
    
    # Adding this line ensures model 1 will always be model 1, 
    # regardless of naming - that way they could be named with numbers
    sf$Model <- factor(sf$Model, levels = c(name1, name2))
    sf$Legend <- factor(sf$Legend, 
                        levels = c(params@species_params$species, "Resource"))
    
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    linecolours <- params@linecolour[legend_levels]
    
    if (power %in% c(0, 1, 2)) {
        y_label <- c("Number density [1/g]", "Biomass density", 
                     "Biomass density [g]")[power + 1]
    }
    else {
        y_label <- paste0("Number density * w^", power)
    }
    
    ggplot(sf, aes(x = w, y = value, colour = Legend, linetype = Model)) +
        geom_line(linewidth = 0.8) +
        scale_x_log10("Weight [g]") +
        scale_y_log10(y_label) + 
        scale_colour_manual(values = linecolours)
    
}

#' Plot the relative difference or percent change between two spectra
#' 
#' This plots a measure of the relative difference between the steady state 
#' spectra of two mizer objects. The user can choose how this difference is 
#' calculated. Let the spectra of the two objects be represented as 
#' \eqn{N_1(w)} and \eqn{N_2(w)}.
#' 
#' If `diff_method` is given as `percent_change`, this function plots the
#' percent change, given by \deqn{ 100*(N_2(w) - N_1(w)) / (N_1(w)).}
#' 
#' If `diff_method` is given as `rel_diff` the difference is calculated 
#' relative to their average, so \deqn{2 (N_2(w) - N_1(w)) / (N_2(w) + N_1(w)).}
#' 
#' The individual spectra are calculated by the [plotSpectra()] function which
#' is passed all additional arguments you supply. So you can for example
#' determine a size range over which to average the simulation results via the
#' `time_range` argument. See [plotSpectra()] for more options.
#' 
#' Note that it does not matter whether the relative difference is calculated
#' for the number density or the biomass density or the biomass density in log
#' weight because the factors of \eqn{w} by which the densities differ cancels 
#' out in the relative difference.
#' 
#' @param object1 An object of class MizerSim or MizerParams
#' @param object2 An object of class MizerSim or MizerParams
#' @param diff_method   The method to calculate the relative change between 
#'                      models. If `percent.change`, the percent change is 
#'                      calculated relative to the value from object 1 with 
#'                      formula 100*(new-old)/old. If `rel.diff` the relative 
#'                      difference is returned given by (new - old)/(old + new).
#' @param power The abundance is plotted as the number density times the weight
#'   raised to this power. The default power = 1 gives the biomass density,
#'   whereas power = 2 gives the biomass density with respect to logarithmic
#'   size bins.
#' @param ... Parameters passed to `plotSpectra()`
#' @concept sumplots
#' @family plotting functions
#' @return A ggplot2 object
#' @export
plotSpectraRelative <- function(object1, object2,
                                power,
                                diff_method = "percent_change", ...) {
    
    sf1 <- mizer::plotSpectra(object1, power = power,
                              return_data = TRUE, ...)
    sf2 <- mizer::plotSpectra(object2, power = power,
                              return_data = TRUE, ...)
    
    # Calculate relative difference
    if (diff_method == "percent_change"){
        sf <-   dplyr::left_join(sf1, sf2, by = c("w", "Legend")) |>
                dplyr::mutate(rel_diff = (value.y - value.x) / value.x)
        yLabel <- "% Change in Biomass"
    } else if (diff_method == "rel_diff"){
        sf <-   dplyr::left_join(sf1, sf2, by = c("w", "Legend")) |>
                dplyr::mutate(rel_diff = ((value.y - value.x) / (value.x + value.y)))
        yLabel <- "Relative Difference in Biomass"
    } else {
        stop("diff_method should be either 'percent_change' or 'rel_diff'.")
    }
    
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    # Ensures that and species stay in the order we expect regardless of
    # naming
    sf$Legend <- factor(sf$Legend, 
                        levels = c(params@species_params$species, "Resource"))
    
    legend_levels <- intersect(names(params@linecolour),
                               unique(sf$Legend))
    linecolours <- params@linecolour[legend_levels]
    
    ggplot(sf, aes(x = w, y = rel_diff, colour = Legend)) +
        geom_line(linewidth = 0.95) +
        labs(x = "Weight [g]", y = yLabel) +
        scale_x_log10() +
        scale_color_manual(values = linecolours) +
        geom_hline(yintercept = 0, linetype = 1,
                   colour = "dark grey", linewidth = 1)
}

#' @rdname plotSpectraRelative
#' @export
plotlySpectraRelative <- function(object1, object2, diff_method,...) {
    ggplotly(plotSpectraRelative(object1, object2, diff_method,...),
             tooltip = c("Legend", "w", "rel_diff"))
}


#' Plot the total productivity for each species Group
#'
#' When called with a \linkS4class{MizerParams} object the total steady 
#' state productivity is plotted for each group. When called with a
#' \linkS4class{MizerSim} object the productivity of each species 
#' through time is plotted.
#' 
#' @inheritSection getProductivity Potential fisheries productivity
#'
#' @param object An object of class \linkS4class{MizerParams} or
#'                  \linkS4class{MizerSim}
#' 
#' @param species   The species to be selected. Optional. By default all
#'                  species are selected. A vector of species names, or a 
#'                  numeric vector with the species indices, or a logical 
#'                  vector indicating for each species whether it is to be 
#'                  selected (TRUE) or not.
#'                  
#' @param end   A boolean value that indicated whether you want the end
#'              productivity of a simulation (default, TRUE) or the
#'              productivity through time. 
#'                  
#' @param total A boolean value that determines whether the total productivity
#'              from all species is plotted as well. Default is FALSE.
#'              
#' @param ylim  A numeric vector of length two providing lower and upper limits
#'              for the y axis. Use NA to refer to the existing minimum or 
#'              maximum. Any values below 1e-20 are always cut off.
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
plotProductivity <- function(object, end  = TRUE,
                             species = NULL, ylim = c(NA, NA),
                             total = FALSE,  return_data = FALSE,
                             min_fishing_l = NULL, max_fishing_l = NULL,
                             start_time = NULL, end_time = NULL,...) {
    
    if (is(object, "MizerSim")) {
        # sim values ----
        sim <- object
        assert_that(is(sim, "MizerSim"),
                    is.flag(return_data))
        
        params <- sim@params
        
        species <- mizer::valid_species_arg(sim, species, 
                                            error_on_empty = TRUE)
        
        if (missing(start_time)) start_time <- 
            as.numeric(dimnames(sim@n)[[1]][1])
        if (missing(end_time)) end_time <- 
            as.numeric(dimnames(sim@n)[[1]][dim(sim@n)[1]])
        if (start_time >= end_time) {
            stop("start_time must be less than end_time")
        }
        
        p <- getProductivity(sim,
                             min_fishing_l = min_fishing_l,
                             max_fishing_l = max_fishing_l, ...)
        # Select time range
        p <- p[(as.numeric(dimnames(p)[[1]]) >= start_time) &
               (as.numeric(dimnames(p)[[1]]) <= end_time), , drop = FALSE]
        
        # Include total
        if (total) {
            p <- cbind(p, Total = rowSums(p))
        }
        
        p <- reshape2::melt(p)
        
        # Implement ylim and a minimal cutoff and bring columns in 
        # desired order
        min_value <- 1e-20
        p <- p[p$value >= min_value &
                     (is.na(ylim[1]) | p$value >= ylim[1]) &
                     (is.na(ylim[2]) | p$value <= ylim[2]), c(1, 3, 2)]
        names(p) <- c("Year", "Biomass", "Species")
        
        if (end == TRUE) { p <- p[end_time, ,drop = TRUE] }
        
        # Select species
        plot_dat <- p[p$Species %in% c("Total", species), ]
        plot_dat$Legend <- plot_dat$Species
        
        if (return_data) return(plot_dat) 
        
        mizer::plotDataFrame(plot_dat, params, xlab = "Year", 
                             ylab = expression(Productivity~"("*g/m^2/year*")"),
                             ytrans = "log10",
                             y_ticks = y_ticks, highlight = highlight,
                             legend_var = "Legend")
        
    } else {
        # params ----
        params <- object
        assert_that(is(params, "MizerParams"),
                    is.flag(return_data))
    
        ### values from object ----
        sp <- params@species_params
        no_sp <- dim(params@interaction)[1]
        
        ### group names ----
        if (is.null(params@species_params$group_names)){
            group_names <- params@species_params$species
            names(group_names) <- params@species_params$species
        } else {
            group_names <- params@species_params$group_names
            names(group_names) <- params@species_params$species
        }
        
        ### get productivity ----
        prod <- getProductivity(params, 
                                min_fishing_l = min_fishing_l,
                                max_fishing_l = max_fishing_l)
        
        ### species selector ----
        sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                    error_on_empty = TRUE)
        species <- dimnames(params@initial_n)$sp[sel_sp]
        species <- gsub('inverts', NA, species)
        species <- species[!is.na(species)]
        sel_sp <- which(!is.na(species))
        prod <- prod[sel_sp, drop = FALSE]
        group_names <- group_names[sel_sp]
        
        ### data frame from selected species ----
        plot_dat <- data.frame(value = prod, Species = species)
        
        ### colors ----
        legend_levels   <- intersect(names(params@linecolour), plot_dat$Species)
        plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
        
        ### return data if requested ----
        if (return_data) return(plot_dat)
        
        # plot ----
        p <- ggplot(plot_dat, aes(x = Species, y = value,
                                  group = Legend, fill = Legend))
        
        p + geom_bar(stat = "identity", position = "dodge") +
            scale_y_continuous(name = expression(Productivity~"("*g/m^2/year*")")) +
            scale_fill_manual(values = params@linecolour[legend_levels],
                              labels = group_names) +
            labs(fill = "Species Group", x = "Species Group")
    }
}


#' @rdname plotProductivity
#' @export
plotlyProductivity <- function(object,
                               species = NULL,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("plot2Productivity", argg),
             tooltip = c("Species", "value"))
}


#' Plot the fisheries productivity of two models or two different size ranges
#'  in the same plot
#'
#' When called with a \linkS4class{MizerParams}
#' object the steady state productivities are plotted.
#'
#' @param object1 First MizerParams or MizerSim object.
#' 
#' @param object2 Second MizerParams or MizerSim object.
#' 
#' @param species   The groups to be selected. Optional. By default all target
#'                  groups are selected. A vector of groups names, or a numeric 
#'                  vector with the groups indices, or a logical vector 
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
                              return_data = FALSE, ...){
    
    # get data frames with plotProductivity ----
    sf1 <- plotProductivity(object1, 
                            species = species,
                            drop = TRUE,
                            min_fishing_l = min_fishing_l1,
                            max_fishing_l = max_fishing_l1,
                            return_data = TRUE, ...)
        sf1$Model <- name1
    sf2 <- plotProductivity(object2, 
                            species = species,
                            drop = TRUE,
                            min_fishing_l = min_fishing_l2,
                            max_fishing_l = max_fishing_l2,
                            return_data = TRUE, ...)
        sf2$Model <- name2
        
    sf <- rbind(sf1, sf2)
    
    # Make sure order of models isnt changed by names
    sf$Model <- factor(sf$Model, levels = c(name1, name2))
    
    # if sim, get params ----
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    # group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }
    
    # plot ----
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    sf$Legend <- factor(sf$Species, levels = legend_levels)
    group_names <- group_names[sf$Legend]
    
    # Return data frame if requested
    if(return_data == TRUE){return(sf)}
    
    if (stack == FALSE){
        p <- ggplot(sf, aes(x = Species, y = value, 
                            group = Model, alpha = Model, fill = Legend))
        
        p + geom_bar(stat = "identity", position = "dodge", color = "black") +
            scale_y_continuous(name = expression(Productivity~"("*g/m^2/year*")")) +
            scale_fill_manual(values = params@linecolour[legend_levels],
                              labels = group_names) +
            scale_alpha_manual(values = c(0.5,1),
                               labels = c(name1, name2)) +
            labs(fill = "Species Group", x = "Species Group") 
        
    } else if (stack == TRUE){
        
        p <- ggplot(sf, aes(x = Model, y = value, 
                            alpha = Model,fill = Legend))
        
        p + geom_bar(stat = "identity", position = "stack", color = "black") +
            scale_y_continuous(name = expression(Productivity~"("*g/m^2/year*")")) +
            scale_fill_manual(values = params@linecolour[legend_levels],
                              labels = group_names) +
            scale_alpha_manual(values = c(0.5,1),
                               labels = c(name1, name2)) +
            labs(fill = "Species Group", x = "Model")
    }
    
}

#' @rdname plot2Productivity
#' @export
plotly2Productivity <- function(object1, object2,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("plot2Productivity", argg),
             tooltip = c("Species", "value"))
}

#' Plot the relative difference between the potential fisheries productivity 
#' rates of two models or two different size ranges in the same plot
#' 
#' This function creates a barplot with the percent change in potential
#' fisheries productivity between either: (1) two different mizerParams objects
#' (2) two different size ranges. If comparing two mizerParams objects, they
#' must have the same species groups. 
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
#' @param species   The groups to be selected. Optional. By default all target
#'                  groups are selected. A vector of groups names, or a numeric 
#'                  vector with the groups indices, or a logical vector 
#'                  indicating for each group whether it is to be selected 
#'                  (TRUE) or not.
#'                  
#' @param diff_method   The method to calculate the relative change between 
#'                      models. If `percent.change`, the percent change is 
#'                      calculated relative to the value from object 1 with 
#'                      formula 100*(new-old)/old. If `rel.diff` the relative 
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
#'          frame with the the productivity for each Species Group
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
                                     return_data = FALSE, ...){
    
    # get data frames with plotProductivity ----
    sf1 <- plotProductivity(object1, 
                            species = species,
                            min_fishing_l = min_fishing_l1,
                            max_fishing_l = max_fishing_l1,
                            return_data = TRUE, ...)
    sf2 <- plotProductivity(object2, 
                            species = species,
                            min_fishing_l = min_fishing_l2,
                            max_fishing_l = max_fishing_l2,
                            return_data = TRUE, ...)
    
    # Calculate relative difference
    if (diff_method == "percent_change"){
        sf <-   dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |> 
                dplyr::mutate(rel_diff = (value.y - value.x) / value.x)
        yLabel <- "% Change in Productivity"
    } else if (diff_method == "rel_diff"){
        sf <-   dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |>
                dplyr::mutate(rel_diff = ((value.y - value.x) / (value.x + value.y)))
        yLabel <- "Relative Difference in Productivity"
    } else {
        stop("diff_method should be either 'percent_change' or 'rel_diff'.")
    }
    
    # return data if requested ----
    if(return_data == TRUE){return(sf)}
    
    # save 1 set for species names ----
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    # group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }

    # plot -----
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    sf$Legend <- factor(sf$Species, levels = legend_levels)
    group_names <- group_names[sf$Legend]
    
    p <- ggplot(sf, aes(x = Species, y = rel_diff,
                        fill = Legend))
    
    p + geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_y_continuous(name = yLabel) +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names) +
        labs(fill = "Species Group", x = "Species Group") +
        geom_hline(yintercept = 0, linetype = 1,
                   colour = "dark grey", linewidth = 0.9)
}

#' @rdname plotProductivityRelative
#' @export
plotlyProductivityRelative <- function(object1, object2, 
                                       diff_method,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("plotProductivityRelative", argg),
             tooltip = c("Species", "value"))
}

#' Plot the total abundance for each species group at steady state
#' 
#' This functions creates a barplot with the abundance of each species group
#' within a given size range. 
#' 
#' @param object An object of class \linkS4class{MizerParams}
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
plotTotalAbundance <- function(object,
                               species = NULL,
                               min_fishing_l = NULL, 
                               max_fishing_l = NULL,
                               return_data = FALSE, ...) {
    
    # object checks ----
    if (is(object, "MizerSim")) {
        ## sim values ----
        # get total abundance at last timestep
        params <- object@params
        end_time  <- max(as.numeric(dimnames(object@n)$time))
        
        abd <- mizer::getN(object,
                           min_l = min_fishing_l,
                           max_l = max_fishing_l)
        
        abd <- abd[end_time, ,drop = TRUE]
        
    } else {
        
    # params ----
    params <- object
    assert_that(is(params, "MizerParams"),
                    is.flag(return_data))
    
    abd <- mizer::getN(params,
                       min_l = min_fishing_l,
                       max_l = max_fishing_l)

    # create plot_dat ----
    ## values from object ----
    sp <- params@species_params
    no_sp <- dim(params@interaction)[1]
    
    ## group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }
    
    ## species selector ----
    sel_sp <- mizer::valid_species_arg(params, species, 
                                       return.logical = TRUE, 
                                       error_on_empty = TRUE)
    species <- dimnames(params@initial_n)$sp[sel_sp]
    species <- species[!is.na(species)]
    sel_sp <- which(!is.na(species))
    abd <- abd[sel_sp, drop = FALSE]
    group_names <- group_names[sel_sp]
    
    ## data frame from selected species ----
    plot_dat <- data.frame(value = abd, Species = species)
    
    ## colors ----
    legend_levels   <- intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
    
    ## return data if requested ----
    if (return_data) return(plot_dat)
    
    # plot ----
    p <- ggplot(plot_dat, aes(x = Species, y = value,
                              group = Legend, fill = Legend))
    
    p + geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(name = expression("Total Abundance (no./m^2)")) +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names) +
        labs(fill = "Species Group", x = "Species Group")
    }
}


#' Plot the total fishable biomass for each Species Group at steady state
#' 
#' This functions creates a barplot with the biomass of each Species Group
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
    
    # object checks ----
    if (is(object, "MizerSim")) {
        ## sim values ----
            # get total biomass at last timestep
        params <- object@params
        end_time  <- max(as.numeric(dimnames(object@n)$time))
        
        biom <- mizer::getBiomass(object, 
                                  min_l = min_fishing_l,
                                  max_l = max_fishing_l)
        
        biom <- biom[end_time, ,drop = TRUE]
    } else {
        
    # params ----
    params <- object
    assert_that(is(params, "MizerParams"),
                is.flag(return_data))
    
    biom <- mizer::getBiomass(params, 
                              min_l = min_fishing_l,
                              max_l = max_fishing_l)
    }
        
    # create plot_dat ----
    ## values from object ----
    sp <- params@species_params
    no_sp <- dim(params@interaction)[1]
    
    ## group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }
    
    ## species selector ----
    sel_sp <- mizer::valid_species_arg(params, species, 
                                       return.logical = TRUE, 
                                       error_on_empty = TRUE)
    species <- dimnames(params@initial_n)$sp[sel_sp]
    species <- species[!is.na(species)]
    sel_sp <- which(!is.na(species))
    biom <- biom[sel_sp, drop = FALSE]
    
    ## data frame from selected species ----
    plot_dat <- data.frame(value = biom, Species = species)
    
    ## colors ----
    legend_levels   <- intersect(names(params@linecolour), plot_dat$Species)
    plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
    
    ## return data if requested ----
    if (return_data){ return(plot_dat) }
    
    # plot ----
    p <- ggplot(plot_dat, aes(x = Species, y = value,
                              group = Legend, fill = Legend))
    
    p + geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(name = expression("Total Biomass"~"("*g/m^2*")")) +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names[legend_levels]) +
        labs(fill = "Species Group", x = "Species Group")
}

#' @rdname plotTotalBiomass
#' @export
plotlyTotalBiomass <- function(object,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("plotTotalBiomass", argg),
             tooltip = c("Species", "value"))
}


#' Plot the total biomass of two models or of two different size ranges in 
#' the same plot
#'
#' When called with a \linkS4class{MizerParams}
#' object the steady state biomasses are plotted.
#' 
#' @inheritParams plot2Productivity
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
plot2TotalBiomass <- function(object1, object2, 
                              species = NULL,
                              name1 = "First", name2 = "Second",
                              min_fishing_l1 = NULL, max_fishing_l1 = NULL,
                              min_fishing_l2 = NULL, max_fishing_l2 = NULL,
                              stack = FALSE,
                              return_data = FALSE, ...){
    
    # get data frames with plotTotalBiomass ----
    sf1 <- plotTotalBiomass(object1, 
                            species = species,
                            min_fishing_l = min_fishing_l1,
                            max_fishing_l = max_fishing_l1,
                            return_data = TRUE, ...)
        sf1$Model <- name1
    sf2 <- plotTotalBiomass(object2, 
                            species = species,
                            min_fishing_l = min_fishing_l2,
                            max_fishing_l = max_fishing_l2,
                            return_data = TRUE, ...)
        sf2$Model <- name2
        
    sf <- rbind(sf1, sf2)
    
    # Make sure model names dont change order
    sf$Model <- factor(sf$Model, levels = c(name1, name2))
    
    # if sim, get params ----
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    # group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }
    
    # plot ----
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    sf$Legend <- factor(sf$Species, levels = legend_levels)
    group_names <- group_names[sf$Legend]
    
    # Return data frame if requested
    if(return_data == TRUE){return(sf)}
    
    if (stack == FALSE){
        p <- ggplot(sf, aes(x = Species, y = value, 
                            group = Model, alpha = Model, fill = Legend))
        
        p + geom_bar(stat = "identity", position = "dodge", color = "black") +
            scale_y_continuous(name = expression("Total Biomass"~"("*g/m^2*")")) +
            scale_fill_manual(values = params@linecolour[legend_levels],
                              labels = group_names) +
            scale_alpha_manual(values = c(0.5,1),
                               labels = c(name1, name2)) +
            labs(fill = "Species Group", x = "Species Group") 
        
    } else if (stack == TRUE){
        
        p <- ggplot(sf, aes(x = Model, y = value, 
                            alpha = Model,fill = Legend))
        
        p + geom_bar(stat = "identity", position = "stack", color = "black") +
            scale_y_continuous(name = expression("Total Biomass"~"("*g/m^2*")")) +
            scale_fill_manual(values = params@linecolour[legend_levels],
                              labels = group_names) +
            scale_alpha_manual(values = c(0.5,1),
                                   labels = c(name1, name2)) +
            labs(fill = "Species Group", x = "Model")
    }
            
}

#' @rdname plot2TotalBiomass
#' @export
plotly2TotalBiomass <- function(object1, object2,
                                species = NULL,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("plot2TotalBiomass", argg),
             tooltip = c("Species", "value"))
}

#' Plot the relative difference in between the total fishable biomasses of each
#' each Species Group at steady state
#' 
#' This functions creates a barplot with the relative change in biomass of 
#' each Species Group within a size range between either (1) two different 
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
                                     return_data = FALSE, ...){
    
    # get data frames with plotTotalBiomass ----
    sf1 <- plotTotalBiomass(object1, 
                            species = species,
                            min_fishing_l = min_fishing_l1,
                            max_fishing_l = max_fishing_l1,
                            return_data = TRUE, ...)
    sf2 <- plotTotalBiomass(object2, 
                            species = species,
                            min_fishing_l = min_fishing_l2,
                            max_fishing_l = max_fishing_l2,
                            return_data = TRUE, ...)
    
    # Calculate relative difference
    if (diff_method == "percent_change"){
        sf <- dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |>
              dplyr::mutate(rel_diff = (value.y - value.x) / value.x)
        
            yLabel <- "% Change in Total Biomass"
            
    } else if (diff_method == "rel_diff"){
        sf <- dplyr::left_join(sf1, sf2, by = c("Species", "Legend")) |>
              dplyr::mutate(rel_diff = ((value.y - value.x) / (value.x + value.y)))
        
            yLabel <- "Relative Difference in Total Biomass"
            
    } else {
        stop("diff_method should be either 'percent_change' or 'rel_diff'.")
    }
    
    # Return data frame if requested
    if(return_data == TRUE){return(sf)}
    
    # if sim, get params ----
    if (is(object1, "MizerSim")) {
        params <- object1@params
    } else {
        params <- object1
    }
    
    # group names ----
    if (is.null(params@species_params$group_names)){
        group_names <- params@species_params$species
        names(group_names) <- params@species_params$species
    } else {
        group_names <- params@species_params$group_names
        names(group_names) <- params@species_params$species
    }
    
    # plot ----
    legend_levels <- intersect(names(params@linecolour), unique(sf$Legend))
    sf$Legend <- factor(sf$Species, levels = legend_levels)
    group_names <- group_names[sf$Legend]
    
    p <- ggplot(sf, aes(x = Species, y = rel_diff, fill = Legend))
    
    p + geom_bar(stat = "identity", position = "dodge", color = "black") +
        scale_y_continuous(name = yLabel) +
        scale_fill_manual(values = params@linecolour[legend_levels],
                          labels = group_names) +
        labs(fill = "Species Group", x = "Species Group") + 
        geom_hline(yintercept = 0, linetype = 1,
                   colour = "dark grey", linewidth = 0.9)
    
}

#' @rdname plotTotalBiomassRelative
#' @export
plotlyTotalBiomassRelative <- function(object1, object2,
                                       diff_method,...) {
    
    argg <- as.list(environment())
    ggplotly(do.call("TotalBiomassRelative", argg),
             tooltip = c("Species", "value"))
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
#' @param object An object of class \linkS4class{MizerParams}
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
#' @import ggplot2
#' @export
#' 
#' @concept sumplots
#' @family plotting functions
#' @seealso [plotTotalAbundance()], [plotTotalBiomass()], [plotProductivity()]
plotRelativeContribution <- function(object,
                                     min_size = NULL,
                                     min_fishing_l = NULL,
                                     return_data = FALSE,...){
    
    abd <- plotTotalAbundance(object,
                              min_fishing_l = min_size,
                              return_data = TRUE, ...)
    abd$Metric <- "Abundance"

    biom <- plotTotalBiomass(object, 
                             min_fishing_l = min_size,
                             return_data = TRUE, ...)
    biom$Metric <- "Biomass"
    
    prod <- plotProductivity(object, 
                             min_fishing_l = min_fishing_l,
                             return_data = TRUE, ...)
    prod$Metric <- "Productivity"
    
    if (is(object, "MizerSim")) { 
        params <- object@params 
    } else {
        params <- object
        assert_that(is(params, "MizerParams"))
    }
    
    # Remove invertebrates
    abd  <- subset(abd,  Species != 'inverts')
    biom <- subset(biom, Species != 'inverts')
    prod <- subset(prod, Species != 'inverts')
    
    # Relative Contribution
    abd  <- dplyr::mutate(abd, rel = value / sum(value) * 100)
    biom <- dplyr::mutate(biom, rel = value / sum(value) * 100)
    prod <- dplyr::mutate(prod, rel = value / sum(value) * 100)
    
    rel <- rbind(abd, biom, prod)

    # Legend       
    legend_levels <- intersect(names(params@linecolour), unique(rel$Legend))
    rel$Legend  <- factor(rel$Species, levels = legend_levels) 
    
    # Return data if requested
    if(return_data == TRUE){ return(rel) }

    # Plot
    p <- ggplot(rel, aes(x = Metric, y = rel, fill = Legend))
    
    p + geom_bar(stat = "identity", position = "fill", color = "black") +
        scale_fill_manual(values = params@linecolour[legend_levels]) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(y = "Relative Contribution", x = "Metric", fill = "")
    
}