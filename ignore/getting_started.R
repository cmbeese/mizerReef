## ----setup, message=FALSE-----------------------------------------------------
library(mizerReef)
data("karpata_species")
data("karpata_int")

# Show the species parameter rows and the three reef-specific columns
if(!all(c('refuge_user','bad_pred','satiation') %in% names(karpata_species))){
  # If the columns are missing, show how to add defaults
  karpata_species$refuge_user <- rep(TRUE, nrow(karpata_species))
  karpata_species$bad_pred   <- rep(FALSE, nrow(karpata_species))
  karpata_species$satiation  <- rep(TRUE, nrow(karpata_species))
}

knitr::kable(karpata_species[ , c('species','refuge_user','bad_pred','satiation')])


## ----refuge-method-bodyshape-examples, fig.height=4, fig.width=10, warning=FALSE----

utils::data("bodyshape_methods_example", package = "mizerReef")
utils::data("bodyshape_methods_example_int", package = "mizerReef")

# Assign to the local variable names used in the example
sp4 <- bodyshape_methods_example
int4 <- bodyshape_methods_example_int

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
                   refuge_density = c(0, 10^5, 10^5, 0, 0,
                                      0, 10^6, 10^3, 0, 0))

## Minimal unstructured-resource interaction for the small example
UR_int <- data.frame(algae = rep(0, nrow(sp4)), detritus = rep(0, nrow(sp4)))

## Create params objects for each method
params_sig <- newReefParams(group_params = sp4,
                            interaction = int4,
                            UR_interaction = UR_int,
                            method = methods[1],
                            method_params = sig,
                            satiation = rep(TRUE, nrow(sp4)),
                            refuge_user = sp4$refuge_user,
                            bad_pred = sp4$bad_pred)
params_sig@other_params$degrade <- FALSE

params_bin <- newReefParams(group_params = sp4,
                            interaction = int4,
                            UR_interaction = UR_int,
                            method = methods[2],
                            method_params = bin,
                            satiation = rep(TRUE, nrow(sp4)),
                            refuge_user = sp4$refuge_user,
                            bad_pred = sp4$bad_pred)
params_bin@other_params$degrade <- FALSE

params_comp <- newReefParams(group_params = sp4,
                             interaction = int4,
                             UR_interaction = UR_int,
                             method = methods[3],
                             method_params = comp,
                             satiation = rep(TRUE, nrow(sp4)),
                             refuge_user = sp4$refuge_user,
                             bad_pred = sp4$bad_pred)
params_comp@other_params$degrade <- FALSE

## Plot with method labels (use plotRefuge and add titles)
g_sig <- plotRefuge(params_sig) + ggplot2::ggtitle('Method: sigmoidal')
g_bin <- plotRefuge(params_bin) + ggplot2::ggtitle('Method: binned')
g_comp <- plotRefuge(params_comp) + ggplot2::ggtitle('Method: competitive')

## Arrange and print
library(ggpubr)
ggpubr::ggarrange(g_sig, g_bin, g_comp, ncol = 3, common.legend = TRUE, legend = 'top')

## Optionally save the figure (example from Chapter 3)
# ggsave('bodyshape_methods_refuge_examples.png', plot = last_plot(), width = 27, height = 12, units = 'cm', dpi = 600)


## ----change-sigmoidal, eval=FALSE---------------------------------------------
# # Update L_refuge and proportion protected for a sigmoidal params object
# newRefuge(params_sig, new_refuge = TRUE, new_L_refuge = 20, new_prop_protect = 0.25)
# # Recalculate refuge and (if needed) retune steady state: reefSteady(), matchReefGrowth(), etc.


## ----change-binned, eval=FALSE------------------------------------------------
# # Scale all bins by 0.5 (50% reduction in protection values)
# newRefuge(params_bin, new_refuge = TRUE, scale_bin = 0.5)
# 
# # For competitive method, scaling refuge_density will change the availability
# # of refuges and therefore the density-dependent allocation.
# newRefuge(params_comp, new_refuge = TRUE, scale_bin = 0.5)


## ----create-run, eval = FALSE-------------------------------------------------
# # Example (eval = FALSE to avoid running a long simulation in CRAN checks)
# params <- newReefParams(group_params = karpata_species,
#                         interaction = karpata_int,
#                         UR_interaction = data.frame(algae = karpata_species$interaction_algae,
#                                                     detritus = karpata_species$interaction_detritus),
#                         method = 'sigmoidal',
#                         method_params = data.frame(L_refuge = 15, prop_protect = 0.2),
#                         refuge_user = karpata_species$refuge_user,
#                         bad_pred = karpata_species$bad_pred,
#                         satiation = karpata_species$satiation)
# 
# # Run a short projection (adjust t_max and effort as needed)
# sim <- project(params, t_max = 5, effort = 1)
# 
# # Explore results
# plot(sim)
# plotRefuge(params)


## ----echo = FALSE-------------------------------------------------------------
sessioninfo::session_info()

