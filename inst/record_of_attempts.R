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

# With these parameters, herbivores consume plankton at small sizes and 
#   transition fully to algae by maturity 
# With these parameters, invertebrates consume plankton and detritus,
#   with the proportion of detritus increasing with size

# Attempt 1 --------------------------------------------------------------------
    
# Won't converge with smaller w_pp_cutoff

    ## Set model
    params <- newReefParams(group_params = karpata_species,
                            interaction = karpata_int,
                            method = "binned", # w_pp_cutoff = 0.1,
                            method_params = tuning_profile)
    
    ## Project to first steady state
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
    # Converges on 7th attempt
    
    ## Calibrate biomasses
    params <- calibrateReefBiomass(params)
    
    # Plot species vs biomass to see which ones are furthest away
    plotBiomassVsSpecies(params)
    
    # Match biomass species by species, starting with ones that are closest
    params <- matchBiomasses(params, species = "parrotfish")
    params <- reefSteady(params)
    
    params <- matchBiomasses(params, species = "farm_damsel")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() 

    params <- matchBiomasses(params, species = "herbs")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    params <- matchBiomasses(params, species = "pred_inv")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    params <- matchBiomasses(params, species = "pred_eng")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    params <- matchBiomasses(params, species = "pred_grab")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_plank")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady()
    
    # Match invertivores again
    params <- matchBiomasses(params, species = "pred_inv")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    # Match biomasses
    params <- calibrateReefBiomass(params)
    params <- matchBiomasses(params)
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady()

    plotBiomassVsSpecies(params)
    
    plotSpectra(params, power = 2)

    params <- setBevertonHolt(params, erepro = 0.0001)   

    # Match observed growth rates
    params <- matchReefGrowth(params)
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()


    # Check for match with age at maturity
    age_mat_observed = karpata_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)
    # Not bad
    
    plotBiomassVsSpecies(params)
    # bad!!
    
    # Iterate to refine biomass
    params <- params |>
        calibrateReefBiomass() |> matchBiomasses()|> reefSteady()|>
        matchReefGrowth()|> reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> reefSteady()|>
        matchReefGrowth()|> reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()
    
    # Error in 
    # mizer::setBevertonHolt(params, 
    #                        reproduction_level = old_reproduction_level) : 
    #     Some species have no reproduction.
    

    rm(params)

# Attempt 2 --------------------------------------------------------------------
# Try with 10 plus biomasses instead of 5 plus - has drastically smaller
# biomass for pred_plank, no farming damselfish
    
    ## Set model
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", # w_pp_cutoff = 0.1,
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
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    # This is still where it gets unhappy
    params <- matchBiomasses(params, species = "pred_plank")
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
    
    rm(params)
    
# Attempt 3 --------------------------------------------------------------------
    # Try with 10 plus biomasses instead of 5 plus - has drastically smaller
    # biomass for pred_plank, no farming damselfish
    # increase beta for planktivores?
    karpata_10plus$beta[karpata_10plus$species == "pred_plank"] <- 1000
    
    ## Set model
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", w_pp_cutoff = 0.1,
                            method_params = tuning_profile)
    
    params <- reefSteady(params)
    
    rm(params)
    
# Attempt 4 --------------------------------------------------------------------
    # Try with 10 plus biomasses instead of 5 plus - has drastically smaller
    # biomass for pred_plank, no farming damselfish
    # increase beta for planktivores?
    karpata_10plus$beta[karpata_10plus$species == "pred_plank"] <- 1000
    karpata_10plus$sigma[karpata_10plus$species == "pred_grab"] <- 2
    
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
        reefSteady() |> reefSteady() |> reefSteady() 
    
    plotBiomassVsSpecies(params)
    
    # Match biomass species by species, starting with ones that are closest
    params <- matchBiomasses(params)
    params <- matchReefGrowth(params)
    params <- reefSteady(params)
    
    params <- reefSteady(params)
    
    # Eventually errors, some species have no reproduction
    
    rm(params)
    
# Attempt 5 --------------------------------------------------------------------
    # Try with 10 plus biomasses instead of 5 plus - has drastically smaller
    # biomass for pred_plank, no farming damselfish
    # increase beta for planktivores?
    karpata_10plus$interaction_algae[karpata_species$interaction_algae == 0.5] <- 1
    karpata_10plus$interaction_algae[karpata_species$interaction_detritus == 0.5] <- 1
    karpata_10plus$beta[karpata_10plus$species == "pred_plank"] <- 1000
    karpata_10plus$sigma[karpata_10plus$species == "pred_grab"] <- 2
    
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
    #params <- matchBiomasses(params)
    #params <- reefSteady(params)
    
    # Error in mizer::setBevertonHolt(params, 
    #                     reproduction_level = old_reproduction_level) : 
    #     Some species have no reproduction.
    
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
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() 
    
    # Error in mizer::setBevertonHolt(params, 
    #                         reproduction_level = old_reproduction_level) : 
    #     Some species have no reproduction.
    
    params <- matchBiomasses(params, species = "pred_inv")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_eng")
    params <- reefSteady(params)
    
    plotBiomassVsSpecies(params)
    
    params <- matchBiomasses(params, species = "pred_grab")
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    plotBiomassVsSpecies(params)
    
    params <- params |>
        matchReefGrowth() |> matchReefGrowth() |> matchReefGrowth() |>
        matchReefGrowth() |> matchReefGrowth()
    
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady()
    
    
    
    
    
    
    
   