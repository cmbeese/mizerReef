# Setting up a Caribbean coral reef mizer model with multiple resources
# Model steady state calibration
# Ten species groups - Karpata reef
# last tuned 26th January 2024

# THIS NO LONGER WORKS
## Setup - load packages -------------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)
library(here)

## Load parameters -------------------------------------------------------------

# Load species parameter data
karpata_species <- read.csv(here("inst/karpata_species.csv"))
karpata_10plus  <- read.csv(here("inst/karpata_10plus.csv"))
karpata_int     <- read.csv(here("inst/cbn_interaction.csv"),  row.names = 1)
karpata_refuge  <- karpata_refuge
tuning_profile  <- tuning_profile

# Attempt 1 --------------------------------------------------------------------
        ## Set model
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", w_pp_cutoff = 0.1,
                            method_params = tuning_profile)
    
    ## Project to first steady state
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
    # Converges on 7th attempt
    
    # Match biomasses
    params <- calibrateReefBiomass(params)
    params <- matchBiomasses(params)
    params <- reefSteady(params)
    
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 17755.7361631237
    # Error in mizer::setBevertonHolt(params, 
    #                                 reproduction_level = old_reproduction_level) :
    #     Some species have no reproduction.
    
    
# Attempt 2 --------------------------------------------------------------------
    
    ## Set model - new set of parameters
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", 
                            method_params = tuning_profile)
    
    ## Project to first steady state
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
    
    # Match biomasses doesn't converge - try species by species
    params <- calibrateReefBiomass(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    # Match biomass species by species, starting with ones that are closest
    params <- matchBiomasses(params, species = "parrotfish")
    params <- reefSteady(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "herbs")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_inv")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_eng")
    params <- reefSteady(params)
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_grab")
    params <- reefSteady(params)
    
    plotBiomassVsSpecies(params)
    
    # This is  where it gets unhappy
    params <- matchBiomasses(params, species = "pred_plank")
    
    # Simulation run did not converge after 99 years. Value returned by 
    # the distance function was: 183067.146045887
    
    
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady()
    
    # Error in mizer::setBevertonHolt(params, 
    #                             reproduction_level = old_reproduction_level) : 
    #     Some species have no reproduction.

    plotPredMort(params)
    plotSpectra(params, power = 2)
    plotDiet(params)
    # Everything is eating planktivores. weird.
    
    rm(params)
    
# Attempt 3 --------------------------------------------------------------------
    # Try with 10 plus biomasses instead of 5 plus - has drastically smaller
    # biomass for pred_plank, no farming damselfish
    # increase beta for planktivores?
    karpata_10plus$beta[karpata_10plus$species == "pred_plank"] <- 1000
    
    ## Set model
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", #w_pp_cutoff = 0.1,
                            method_params = tuning_profile)
    
    ## Project to first steady state
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
    
    # Match biomasses doesn't converge - try species by species
    params <- calibrateReefBiomass(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    # Match biomass species by species, starting with ones that are closest
    params <- matchBiomasses(params, species = "parrotfish")
    params <- reefSteady(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "herbs")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_inv")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_eng")
    params <- reefSteady(params)
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_grab")
    params <- reefSteady(params)
    
    plotBiomassVsSpecies(params)
    
    # This is  where it gets unhappy
    params <- matchBiomasses(params, species = "pred_plank")
    params <- reefSteady(params)
    
    plotPredMort(params)
    plotSpectra(params, power = 2)
    plotDiet(params)
    # Looking slightly better, lets increase again 
    rm(params)
    
# Attempt 4 --------------------------------------------------------------------
    # Try with 10 plus biomasses instead of 5 plus - has drastically smaller
    # biomass for pred_plank, no farming damselfish
    # increase beta for planktivores?
    karpata_10plus$beta[karpata_10plus$species == "pred_plank"] <- 8000
    
    ## Set model
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", #w_pp_cutoff = 0.1,
                            method_params = tuning_profile)
    
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
    
    ## Calibrate biomasses
    params <- calibrateReefBiomass(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    # Match biomass species by species, starting with ones that are closest
    params <- matchBiomasses(params, species = "parrotfish")
    params <- reefSteady(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "herbs")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_eng")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() 
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_grab")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_inv")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    # This is  where it gets unhappy
    params <- matchBiomasses(params, species = "pred_plank")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 506032.221535751
    # Error in mizer::setBevertonHolt(params, 
    #                                 reproduction_level = old_reproduction_level) : 
    #     Some species have no reproduction.
    
    plotPredMort(params)
    plotSpectra(params, power = 2)
    plotDiet(params)
    # Looking slightly better, lets increase again 
    rm(params)
    
# Attempt 5 --------------------------------------------------------------------
    # Reduce interactions with planktivores?
    karpata_10plus$beta[karpata_10plus$species == "pred_plank"] <- 10000
    karpata_int$pred_plank[karpata_int$pred_plank == 1] <- 0.5
    
    ## Set model
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", #w_pp_cutoff = 0.1,
                            method_params = tuning_profile)
    
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
    
    ## Calibrate biomasses
    params <- calibrateReefBiomass(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    # Match biomass species by species, starting with ones that are closest
    #params <- matchBiomasses(params)
    #params <- reefSteady(params)
    
    
    # Match biomass species by species, starting with ones that are closest
    params <- matchBiomasses(params, species = "parrotfish")
    params <- reefSteady(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "herbs")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() 
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_plank")
    params <- reefSteady(params)
    
    params <- matchBiomasses(params, species = "pred_inv")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_eng")
    params <- reefSteady(params)
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_grab")
    params <- reefSteady(params) 
    plotBiomassVsSpecies(params)
    
    # Now all together
    params <- calibrateReefBiomass(params)
    params <- matchBiomasses(params)
    params <- reefSteady(params)
    
    # Now match growth
    params <- matchReefGrowth(params)
    
    # With these parameter values the pred_grab does not have enough food 
    # to cover its metabolic cost NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
    plotPredMort(params)
    plotSpectra(params, power = 2)
    plotDiet(params)
    # Looking slightly better, lets increase again 
    rm(params)
    
# Attempt 6 --------------------------------------------------------------------
    # Reduce interactions with planktivores?
    karpata_10plus$beta[karpata_10plus$species == "pred_plank"] <- 10000
    karpata_int$pred_plank[karpata_int$pred_plank == 1] <- 0.5
    
    ## Set model
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", w_pp_cutoff = 0.1,
                            method_params = tuning_profile)
    
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
    
    ## Calibrate biomasses
    params <- calibrateReefBiomass(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    # Match biomass species by species, starting with ones that are closest
    # params <- matchBiomasses(params)
    # params <- reefSteady(params)
    # 
    
    # Match biomass species by species, starting with ones that are closest
    params <- matchBiomasses(params, species = "parrotfish")
    params <- reefSteady(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "herbs")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() 
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_plank")
    params <- reefSteady(params)
    
    params <- matchBiomasses(params, species = "pred_inv")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    # OSCILLATING!
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 6.84237482032618
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 1430.33280917617
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 2278.06090068667
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 8.08335509847775
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 1737.67642915773
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 1505.7909620901
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 7043.75462585284
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_eng")
    params <- reefSteady(params)
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_grab")
    params <- reefSteady(params) 
    plotBiomassVsSpecies(params)
    
    # Now all together
    params <- calibrateReefBiomass(params)
    params <- matchBiomasses(params)
    params <- reefSteady(params)
    
    # Now match growth
    params <- matchReefGrowth(params)
    
    # With these parameter values the pred_grab does not have enough food 
    # to cover its metabolic cost NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    
    plotPredMort(params)
    plotSpectra(params, power = 2)
    plotDiet(params)
    # Looking slightly better, lets increase again 
    rm(params)
    
    
    
    
    
   