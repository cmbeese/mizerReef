# Setting up a Caribbean coral reef mizer model with multiple resources
# Three groups: Predators, Herbivores, Invertebrates
# Model steady state calibration
# last tuned 20th February 2024

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
constant_tune   <- constant_tune
step_tune       <- step_tune

# Increase refuge in tuning profile
#step_tune$prop_protect <- 2*step_tune$prop_protect
step_tune$prop_protect <- 0.5*step_tune$prop_protect
# With these parameters, herbivores consume plankton at small sizes and 
#   transition fully to algae by maturity 
# With these parameters, invertebrates consume plankton and detritus,
#   with the proportion of detritus increasing with size

## Set model -------------------------------------------------------------------
    params <- newReefParams(group_params = bonaire_species,
                            interaction = bonaire_int,
                            method = "binned",
                            method_params = step_tune)
                            # method_params = tuning_profile)

## Project to first steady state -----------------------------------------------
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |>reefSteady() 

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
    
    # Check biomass match - still way off
    plotBiomassVsSpecies(params)

    # Iterate to refine biomass - run this three times
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
    
## Now switch to competitive method --------------------------------------------
    params <- newRefuge(params,
                        new_method = "competitive",
                        new_method_params = bonaire_refuge)
    
    # Match biomasses again
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
    
    plotBiomassVsSpecies(params) # spot on
    
    # Check match with observed age at maturity
    age_mat_observed = bonaire_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)
    # Also spot on

## Check resulting spectra and tune resources ----------------------------------

    # Resource looks too high - should match sheldon's spectrum
    # looks fairly straight not bad but some bumps
    plotSpectra(params, total = TRUE, power = 1)
    plotSpectra(params, total = TRUE, power = 2)
    
    # plot feeding level to check if resource is too low
    plotFeedingLevel(params, species = "inverts")
    
    # Invertebrate feeding level is stable throughout life - there is enough
    # reosurce
    
# Tune reproduction ------------------------------------------------------------
    # We do not have yield or catch data - can't tune size distribution
    # First attempt to set very low to see what the minimum values are
    params <- setBevertonHolt(params, erepro = 0.0001)
    # Now set setting erepro same for all species, as low as possible
    params <- setBevertonHolt(params, erepro = 0.06)
    params <- reefSteady(params)
    # Check reproduction level (value between 0 and 1) - should be higher for
    # larger, slow growing species and low for small, fast growing ones
    rep <- getReproductionLevel(params)
    # These look good. A reproduction level closer to one means reproduction 
    # rate is almost totally independent of the investment into reproduction
    # Reproduction should be density independent on reefs
    
    # Check comparison of density dependent & independent reduction
    getRDI(params) / getRDD(params)
    # Reproduction is equally density independent nad density dependent for 
    # invwerts, more density independent for herbivores
    
    # Let's increase reproduction level to 0.5 for predators and herbivores
    rep_level <- c(0.5, 0.5, 0.5)
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
    
    # Check new spectra
    plotSpectra(params, total = TRUE, power = 1)
    plotSpectra(params, total = TRUE, power = 2)

# Plots ------------------------------------------------------------------------
    plotBiomassVsSpecies(params)
    plotRefuge(params)
    plotSpectra(params, power = 1, total = TRUE)
    plotDiet(params)  
    plotGrowthCurves(params)
    plotPredMort(params)

    # Save!
    bon_test4 <- reefSteady(params)
    bon_species4 <- bonaire_species

# Save in package --------------------------------------------------------------
    # Params object
    save(bon_test4,   file = "data/bon_test4.rda")
    
    # CSV Files
    save(bon_species4,    file = "data/bon_species4.rda")
    save(bonaire_int,     file = "data/bonaire_int.rda")
    
    # Things that dont change
    save(bonaire_refuge, file = "data/bonaire_refuge.rda")
    save(constant_tune,  file = "data/constant_tune.rda")
    save(step_tune,      file = "data/step_tune.rda")
    