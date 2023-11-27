# Setting up a Caribbean coral reef mizer model with multiple resources

#### Setup - loading packages and functions ------------------------------------
# Load in relevant packages
library(mizer)
library(mizerReef)
library(assertthat)

#### Parameters ----------------------------------------------------------------
# Load species parameter data
cbn_species <- read.csv("data/cbn_species.csv", row.names = 1, header = TRUE)
cbn_int     <- read.csv("data/cbn_int.csv",     row.names = 1, header = TRUE)
cbn_UR_int  <- read.csv("data/cbn_UR_int.csv",  row.names = 1, header = TRUE)

# Create some tester refuge scenarios
methods <- c("sigmoidal", "binned", "competitive")
rogers_2014 <- data.frame(L_refuge = 15, prop_protect = 0.2)
beese_2023  <- data.frame(start_L = seq(0, 45, 5),
                          end_L = seq(5, 50, 5),
                          prop_protect = c(1, 1, 0.9, rep(0,7)))
rogers_2018 <- data.frame(start_L = seq(0, 45, 5),
                          end_L = seq(5, 50, 5),
                          refuge_density = c(0, 10^5, 10^5, 0, 0))

## SET MODEL -------------------------------------------------------------------
cbn_1 <- newReefParams(species_params = cbn_species,
                       interaction = cbn_int,
                       UR_interaction = cbn_UR_int,
                       method = methods[1],
                       method_params = rogers_2014)

cbn_2 <- newReefParams(species_params = cbn_species,
                            interaction = cbn_int,
                            UR_interaction = cbn_UR_int,
                            method = methods[2],
                            method_params = beese_2023)

cbn_3 <- newReefParams(species_params = cbn_species,
                            interaction = cbn_int,
                            UR_interaction = cbn_UR_int,
                            method = methods[3],
                            method_params = rogers_2018)

# Refuge plots
one <- plotRefuge(cbn_1)
two <- plotRefuge(cbn_3)
thr <- plotRefuge(cbn_2)

print(one)
print(two)
print(thr)

# Save as rda
save(cbn_params, file = "data/cbn_params.rda")

## Project to steady state
cbn_params <- reef_steady(cbn_params)

# Plots
plotVulnerable(cbn_params)
plotSpectra(cbn_params, power = 2)
plotBiomass(cbn_params)

# ## FOR ERRORS ------------------------------------------------------------------
# rates_fns <- lapply(cbn_sim@rates_funcs, get)
# 
# r <- rates_fns$Rates(cbn_sim,
#                      n = cbn_sim@initial_n,
#                      n_pp = cbn_sim@initial_n_pp,
#                      n_other = cbn_sim@initial_n_other,
#                      t = 0,
#                      effort = cbn_sim@initial_effort,
#                      rates_fns = rates_fns)
# 
# g <- newMultispeciesParams(NS_species_params)
