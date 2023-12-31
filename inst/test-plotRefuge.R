plotRefuge <- function(object, 
                       species = NULL,
                       all.sizes = FALSE,
                       return_data = FALSE,
                       ...) {
    
    assert_that(is.flag(all.sizes),
                is.flag(return_data))
    
    if (is(object, "MizerSim")) {
        
        stop('This functionality is not set up yet you dumbass.')
        
    } else if (is(object, "MizerParams")) {
        
        params <- object
        sp <- params@species_params
        no_sp <- dim(params@interaction)[1]
        
        if (is.null(sp$group_names)){
            group_names <- sp$species
        } else {
            group_names <- sp$group_names
        }
        
        # Calculate proportion of fish in refuge
        vul <- getVulnerable(params)
        refuge <- (1-vul)
        
        # selector for desired species
        sel_sp <- valid_species_arg(params, species, return.logical = TRUE,
                                    error_on_empty = TRUE)
        species <- dimnames(params@initial_n)$sp[sel_sp]
        species <- gsub('inverts',NA,species)
        species <- species[!is.na(species)]
        species <- 
        refuge <- refuge[species, ,drop = FALSE]
        
        group_length_bins <- matrix(0, nrow = no_sp, ncol = length(params@w))
        for (i in 1:no_sp) {
            group_length_bins[i,] <- (params@w / sp$a[i])^(1 / sp$b[i])
        }
        group_length_bins <- group_length_bins[sel_sp, , drop = FALSE]
        
        max_length = max(sp$l_max)
        
        # Make data from from selected species
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
        
        if (return_data) return(plot_dat)
        
        # Set up colors
        legend_levels <- intersect(names(params@linecolour), plot_dat$Species)
        plot_dat$Legend <- factor(plot_dat$Species, levels = legend_levels)
        plot_dat$Species <- factor(plot_dat$Species, levels = legend_levels)
        linesize <- rep(0.8, length(legend_levels))
        names(linesize) <- names(params@linetype[legend_levels])
        
        p <- ggplot(plot_dat, aes(group = Species)) +
            facet_wrap(~ Species, scales = "free_x") +
            theme(strip.text.x = element_text(size = 6))
        #    strip.background = element_blank(),
        #   strip.text.x =element_blank())
        
        p + geom_line(aes(x = l, y = value,
                          colour = Legend, linetype = Legend,
                          linewidth = Legend)) +
            labs(colour = 'Functional Group', linetype = 'Functional Group',
                 linewidth = 'Functional Group') +
            scale_x_continuous(name = "Total Length [cm]",
                               limits = c(0,max_length)) +
            #scale_x_continuous(name = "Log Size [g]", trans = "log10",
                               #breaks = c(10^-2, 10^0, 10^2, 10^4),
                               #labels = c(-2, 0, 2, 4)) +
            scale_y_continuous(name = "Proportion Protected", limits = c(0, 1)) +
            scale_colour_manual(values = params@linecolour[legend_levels],
                                labels = group_names) +
            scale_linetype_manual(values = params@linetype[legend_levels],
                                  labels = group_names) +
            scale_discrete_manual("linewidth", values = linesize,
                                  labels = group_names)
    }
}


plotRefuge(com_params)
