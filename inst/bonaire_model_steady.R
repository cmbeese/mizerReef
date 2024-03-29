# Setting up a Caribbean coral reef mizer model with multiple resources
# Three groups: Predators, Herbivores, Invertebrates
# Model steady state calibration
# last tuned 23rd February 2024

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
bonaire_refuge  <- bonaire_refuge # Average refuge density per square meter 
tuning_profile  <- tuning_profile # 60% refuge for all size classes

# With these parameters, herbivores consume plankton at small sizes and 
#   transition fully to algae by maturity 
# With these parameters, invertebrates consume plankton and detritus,
#   with the proportion of detritus increasing with size

## Set model -------------------------------------------------------------------
    params <- newReefParams(group_params = bonaire_species,
                            interaction = bonaire_int,
                            method = "binned",
                            method_params = tuning_profile)

## Project to first steady state -----------------------------------------------
    params <- reefSteady(params)

## Calibrate biomasses and growth ----------------------------------------------

    # Match observed species group biomasses
    params <- calibrateReefBiomass(params)
    params <- matchBiomasses(params)
    params <- reefSteady(params)
    
    # Match observed growth rates
    params <- matchReefGrowth(params)
    params <- reefSteady(params)
    
    # Check for match with age at maturity
    age_mat_observed = bonaire_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)
    # Not bad
    
    # Check biomass match
    plotBiomassVsSpecies(params)
    # Biomasses way off from observations
    
    # Iterate to refine biomass - run twice
    params <- params |>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()
    
    plotBiomassVsSpecies(params) # spot on

    # Check match with observed age at maturity
    age_mat_observed = bonaire_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)
    # Closer than needed
    
    plotTotalAbundance(params)
    plotTotalBiomass(params)
    
## Now switch to competitive method --------------------------------------------
    params <- newRefuge(params,
                        new_method = "competitive",
                        new_method_params = bonaire_refuge)
    
    # Match biomasses again - run three times
    params <- params |>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()
    
    # Make sure new refuge is in place
    plotVulnerable(params)
    plotRefuge(params)
    # looks good
    
    plotBiomassVsSpecies(params) # spot on
    
    # Check match with observed age at maturity
    age_mat_observed = bonaire_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)
    # Also spot on

## Check resulting spectra and tune resources ----------------------------------
    
    # Spectra should be reasonably straight to match predictions of Sheldon's
    # spectrum but also have nonlinearities at refuge sizes
    plotSpectra(params, total = TRUE, power = 1) # looks straight, some bumps
    plotSpectra(params, total = TRUE, power = 2) # resource looks low
    
    # plot feeding level to check if resource is too low
    plotFeedingLevel(params, species = "inverts")
    # Invertebrate feeding level is stable throughout life - there is enough
    # resource to support fish, not too little or too much
    
# Tune reproduction ------------------------------------------------------------
    # We do not have yield or catch data - can't tune size distribution
    # First attempt to set very low to see what the minimum values are
    params <- setBevertonHolt(params, erepro = 0.0001)
    # Now set setting erepro same for all species, as low as possible
    params <- setBevertonHolt(params, erepro = 0.009)
    params <- reefSteady(params)
    # Check reproduction level (value between 0 and 1) - should be higher for
    # larger, slow growing species and low for small, fast growing ones
    rep <- getReproductionLevel(params)
    # These are very low for predators, should be higher and too high
    # for invertebrates. A reproduction level closer to one means reproduction 
    # rate is almost totally independent of the investment into reproduction
    # Reproduction should be somewhat density independent on reefs
    
    # Check comparison of density dependent & independent reduction
    getRDI(params) / getRDD(params)
    # preds, 1:1 - maybe too density dependent 
    # inverts, herbs 20+:1 more density independent - maybe too much
    
    # increase reproduction level to 0.5 for all
    rep_level <- c(0.5, 0.5, rep["inverts"])
    names(rep_level) <- c("predators","herbivores","inverts")
    params <- setBevertonHolt(params,
                              reproduction_level = rep_level)
    
    # Iterate to get back to steady state
    params <- params |>
        reefSteady()|>
        reefSteady()|>
        reefSteady()
    
    # Check new reproduction - these look better
    rep <- getReproductionLevel(params)
    getRDI(params) / getRDD(params)
    # Now density independent reproduction is double for predators
    # And reproduction is somewhat density dependent for herbs and inverts
    
    # Check new spectra
    plotSpectra(params, total = TRUE, power = 1)
    plotSpectra(params, total = TRUE, power = 2)
    # These still look good

# Plots ------------------------------------------------------------------------
    plotTotalAbundance(params) # Total abundances look reasonable, inverts in range
    plotTotalBiomass(params)
    plotBiomassVsSpecies(params)
    plotRefuge(params)
    plotSpectra(params, power = 1, total = TRUE)
    plotDiet(params)  
    plotGrowthCurves(params)
    plotPredMort(params)
    # Everything looks good here! I am happy with my results.
    
    # Save!
    bonaire_model <- reefSteady(params)

# Save in package --------------------------------------------------------------
    # 
    # bonaire_model <- bonaire_model
    # bonaire_model@other_params$degrade <- FALSE
    # bonaire_model@other_params$new_refuge <- FALSE
    
    # Params object
    save(bonaire_model,   file = "data/bonaire_model.rda")
    
    # CSV Files
    save(bonaire_species,    file = "data/bonaire_species.rda")
    save(bonaire_int,        file = "data/bonaire_int.rda")
    
    # Things that dont change
    save(bonaire_refuge,  file = "data/bonaire_refuge.rda")
    save(tuning_profile,  file = "data/tuning_profile.rda")
    