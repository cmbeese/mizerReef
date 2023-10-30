# Setting up a Caribbean coral reef mizer model with multiple resources

#### Setup - loading packages and functions ------------------------------------
# Load in relevant packages
library(mizerReef)
library(assertthat)

#### Parameters ----------------------------------------------------------------
# Load species parameter data
rogers_species <- read.csv("data/rogers_species.csv", header = TRUE)
rogers_int     <- read.csv("data/rogers_int.csv",     row.names = 1, header = TRUE)
rogers_UR_int  <- read.csv("data/rogers_UR_int.csv",  row.names = 1, header = TRUE)

# Create some tester refuge scenarios
met <- c("simple", "binned", "data")
rogers_2014 <- data.frame(max_L = 15, prop_protect = 0.2)

## SET MODEL -------------------------------------------------------------------
rogers_params <- newReefParams(species_params = rogers_species,
                               interaction = rogers_int,
                               UR_interaction = rogers_UR_int,
                               method = met[1],
                               method_params = rogers_2014)

# Save as rda
save(rogers_params, file = "data/rogers_params.rda")

## Project to steady state
rogers_params <- reef_steady(rogers_params)

# Plots ------------------------------------------------------------------------
plotVulnerable(rogers_params)
plotSpectra(rogers_params, power = 2)
plotBiomass(rogers_params)

