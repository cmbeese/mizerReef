# Assign to the local variable names used in the example
sp  <- body_shape_example_species_params
int <- body_shape_example_int

## Methods and method parameters used in Chapter 3
methods <- c("sigmoidal", "binned", "competitive")

## Sigmoidal example (one-line specification)
sig <- data.frame(L_refuge = 15, prop_protect = 0.2)

## Binned example: start_L, end_L (cm) with prop_protect per bin
bin <- data.frame(start_L = seq(0, 45, 5),
                  end_L   = seq(5, 50, 5),
                  prop_protect = c(1.0, 0.8, 0.6, 0.4, 0.4,
                                   0.4, 0.3, 0.3, 0, 0))

## Competitive example: refuge density (no./m^2) per bin
comp <- data.frame(start_L = seq(0, 45, 5),
                   end_L   = seq(5, 50, 5),
                   refuge_density = c(190, 10^6, 10^7, 10^4, 0,
                                      10^4, 10^4.5, 100, 1000, 0))

## Minimal unstructured-resource interaction for the small example
UR_int <- data.frame(algae = rep(0, nrow(sp)),
                     detritus = rep(0, nrow(sp)))

## Create params objects for each method
params_sig <- newReefParams(group_params = sp,
                            interaction = int,
                            UR_interaction = UR_int,
                            method = methods[1],
                            method_params = sig,
                            satiation = rep(TRUE, nrow(sp)),
                            refuge_user = sp$refuge_user,
                            bad_pred = sp$bad_pred)
params_sig@other_params$degrade <- FALSE

params_bin <- newReefParams(group_params = sp,
                            interaction = int,
                            UR_interaction = UR_int,
                            method = methods[2],
                            method_params = bin,
                            satiation = rep(TRUE, nrow(sp)),
                            refuge_user = sp$refuge_user,
                            bad_pred = sp$bad_pred)
params_bin@other_params$degrade <- FALSE

params_comp <- newReefParams(group_params = sp,
                             interaction = int,
                             UR_interaction = UR_int,
                             method = methods[3],
                             method_params = comp,
                             satiation = rep(TRUE, nrow(sp)),
                             refuge_user = sp$refuge_user,
                             bad_pred = sp$bad_pred)
params_comp@other_params$degrade <- FALSE
g_sig <- plotRefuge(params_sig) +
    ggplot2::ggtitle('Sigmoidal')
g_bin <- plotRefuge(params_bin) +
    ggplot2::ggtitle('Binned')
g_comp <- plotRefuge(params_comp) +
    ggplot2::ggtitle('Competitive')

print(g_comp)
## Arrange and print
library(ggpubr)
ggpubr::ggarrange(g_sig, g_bin, g_comp,
                  ncol = 3, common.legend = TRUE, legend = 'top')

ggsave('bodyshape_methods_refuge_examples.png',
       plot = last_plot(), width = 27, height = 12,
      units = 'cm', dpi = 600)
