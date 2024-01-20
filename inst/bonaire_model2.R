# Setting up a Caribbean coral reef mizer model with multiple resources
# Model steady state calibration
# last tuned 21st January 2024

## Setup - load packages -------------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)
library(here)

## Load parameters -------------------------------------------------------------

# Load species parameter data
bonaire_species <- bonaire_species
bonaire_int     <- bonaire_int
bonaire_refuge  <- bonaire_refuge

# With these parameters, herbivores only consume algae
# With these parameters, invertebrates only consume detritus
bonaire_species$interaction_resource <- c(1,0,1)

## Set model -------------------------------------------------------------------
params <- newReefParams(group_params = bonaire_species,
                        interaction = bonaire_int,
                        method = "competitive",
                        method_params = bonaire_refuge)

# Fix gammas - this will be internal but not for this model
params@species_params$gamma <- c(6.4,0.2,0.2)

## Project to steady state -----------------------------------------------------
params <- params |>
    reef_steady() |> reef_steady() |> 
    reef_steady() |> reef_steady() |>
    reef_steady() |> reef_steady() 

## Calibrate biomasses and growth ----------------------------------------------

    # Match observed species group biomasses
    params <- calibrateReefBiomass(params)
    params <- matchBiomasses(params)
    params <- reef_steady(params)
    
    # Match observed growth rates
    params <- matchReefGrowth(params)
    params <- reef_steady(params)
    
    # Check for match with age at maturity
    age_mat_observed = bonaire_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)

## Steady state iteration ------------------------------------------------------

    # Iterate to refine growth and biomass
    params <- params |>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()
    
    # Check match with observed age at maturity
    age_mat_observed = bonaire_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)

## Check resulting spectra and tune resources ----------------------------------

# Resource looks low - should match sheldon's spectrum
# looks fairly straight not bad but some bumps
plotSpectra(params, total = TRUE, power = 1)
plotSpectra(params, total = TRUE, power = 2)

# plot feeding level to check if resource is too low
plotFeedingLevel(params, species = "inverts")

# Invert feeding level is relatively stable through life, non-linearities
#   are probably due to refuge

# Tune reproduction ------------------------------------------------------------
# Using resilience to fishing
params <- setBevertonHolt(params, erepro = 0.0001)

params <- reef_steady(params)

# and plot F curves again
plotYieldVsF(params, species = "Predators", 
             F_range = seq(0.1, 0.5, 0.02))

# Cephalopholis cruentata medium resilience species, should have MSY in
# range 0.1 to 0.5
plotYieldVsF(cel_model, species = "Predators", 
             F_range = seq(0.1, 0.9, 0.02))

# Sparisoma viride medium resilience species, should have MSY in
# range 0.1 to 0.5
plotYieldVsF(cel_model, species = "Herbivores", 
             F_range = seq(0.1, 0.9, 0.02))

# Plots ------------------------------------------------------------------------
plotBiomassVsSpecies(params)
plotRefuge(params)
plotSpectra(params, power = 1, total = TRUE)

# Save object ------------------------------------------------------------------

# Params object
save(params,   file = "data/bonaire_model.rda")

# CSV Files
save(bonaire_species, file = "data/bonaire_species.rda")
save(bonaire_int,     file = "data/bonaire_int.rda")
save(bonaire_refuge,  file = "data/bonaire_refuge.rda")

# Build website ----------------------------------------------------------------
pkgdown::build_site()

