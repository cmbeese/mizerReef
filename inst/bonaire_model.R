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
    bonaire_species <- read.csv(here("inst/bonaire_species.csv"))
    bonaire_int     <- read.csv(here("inst/bonaire_int.csv"),  row.names = 1)
    bonaire_refuge  <- bonaire_refuge
    
    # With these parameters, herbivores consume plankton at small sizes and 
    #   transition fully to algae by maturity 
    # With these parameters, invertebrates consume plankton and detritus,
    #   with the proportion of detritus increasing with size

## Set model -------------------------------------------------------------------
    params <- newReefParams(group_params = bonaire_species,
                            interaction = bonaire_int,
                            method = "competitive",
                            method_params = bonaire_refuge)
    
    # Fix gammas - this will be internal but not for this model
    params@species_params$gamma <- c(6.4,0.2,0.2)
    
## Project to steady state -----------------------------------------------------
    params <- params |>
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
    # We do not have yield or catch data - can't tune size distribution
    # Using resilience to fishing following Ken Andersen 2016 procedure
    # Start by setting erepro same for all species
    params <- setBevertonHolt(params, erepro = 0.25)
    # Check reproduction level (value between 0 and 1) - should be higher for
    # larger, slow growing species and low for small, fast growing ones
    getReproductionLevel(params)
    # These look good
    
    # Check comparison of density dependent & independent reduction
    getRDI(params) / getRDD(params)
    # Reproduction is density independent for inverts, density dependent for 
    # preds and herbs
    
    # Iterate to get back to steady state
    params <- params |>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()

    params <- reef_steady(params)
# Plots ------------------------------------------------------------------------
plotBiomassVsSpecies(params)
plotRefuge(params)
plotSpectra(params, power = 1, total = TRUE)
plotDiet(params)  

# I am happy with these parameters!
# Save object ------------------------------------------------------------------

    # Params object
    save(params,   file = "data/bonaire_model.rda")

    # CSV Files
    save(bonaire_species, file = "data/bonaire_species.rda")
    save(bonaire_int,     file = "data/bonaire_int.rda")
    save(bonaire_refuge,  file = "data/bonaire_refuge.rda")

