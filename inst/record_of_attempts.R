# Setting up a Caribbean coral reef mizer model with multiple resources
# Model steady state calibration
# Ten species groups - Karpata reef
# last tuned 26th January 2024

# THIS NO LONGER WORKS
## Setup - load packages -------------------------------------------------------
library(ggplot2)
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)
library(here)

## Load parameters -------------------------------------------------------------

rm(params, karpata_10plus, karpata_int)

# Load species parameter data
# karpata_species <- read.csv(here("inst/karpata_species.csv"))
karpata_10plus  <- read.csv(here("inst/karpata_10plus.csv"))
karpata_int     <- read.csv(here("inst/cbn_interaction.csv"),  row.names = 1)
karpata_refuge  <- karpata_refuge
tuning_profile  <- tuning_profile

# Attempt 1 --------------------------------------------------------------------
    # Remove growth parameters that mizer struggles to match
    #karpata_10plus$age_mat  <- NULL
    karpata_10plus$k_vb     <- NULL
    karpata_10plus$l_mat    <- NULL
    karpata_10plus$interaction_resource[karpata_10plus$species == "pred_inv"] <- 0.3
    #karpata_10plus$beta[karpata_10plus$species == "pred_inv"] <- 30
    karpata_10plus$sigma[karpata_10plus$species == "pred_eng"] <- 1.5

    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", w_pp_cutoff = 0.9,
                            crit_feed = 0.85,
                            method_params = tuning_profile)
    
    rdi <- rep(0.8, dim(karpata_int)[1])
    
    params <- setBevertonHolt(params, reproduction_level = rdi)
    
    ## Project to first steady state
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() 
    # Converges on 7th attempt
    
    # Match biomasses
    params <- calibrateReefBiomass(params)
    params <- matchBiomasses(params)
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() 
    
    plotBiomassVsSpecies(params)
    # Spot on
    
    # Now match growth
    params <- matchReefGrowth(params)
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
    
    
    plotSpectra(params, power = 2, total = TRUE)
    plotFeedingLevel(params)
    
    # # Scale down background??
    params <- scaleReefBackground(params, 2)
    
    params <- params |>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady() |> reefSteady()

    
    # Check for match with age at maturity
    age_mat_observed = karpata_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)
    # Not great
    
    plotBiomassVsSpecies(params)
    # Spot on
    
    params <- setBevertonHolt(params, erepro = 0.0001)
    getReproductionLevel(params)
    
    params <- setBevertonHolt(params, erepro = 0.008)
    getReproductionLevel(params)
    
    #params <- setBevertonHolt(params, species =  reproduction_level = rdi)
    
    params <- reefSteady(params)
    
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady()
    
    getReproductionLevel(params)
    
    
    plotSpectra(params, power = 2)
    plotPredMort(params)
    plotFeedingLevel(params)
    plotDiet(params) + scale_x_log10(limits = c(1, 10000))
    plotDiet(karpata_model2)+ scale_x_log10(limits = c(1, 10000))
    params <- tuneParams(params)
    
    ## Now switch to competitive method --------------------------------------------
    params <- newRefuge(params,
                        new_method = "competitive",
                        new_method_params = karpata_refuge)
    
    # Match biomasses again
    params <- params |>
        matchBiomasses()|> reefSteady()|> 
        matchBiomasses()|> reefSteady()|>
        matchBiomasses()|> reefSteady()|>
        matchBiomasses()|> reefSteady()|> reefSteady()|>
        reefSteady()|> reefSteady()|> reefSteady()|> reefSteady()
    
    # Make sure new refuge is in place
    plotRefuge(params)
    
    plotBiomassVsSpecies(params) # spot on
    
    plotSpectra(params, power = 2, total = TRUE)
    
    
    # Save!
    karpata_model <- reefSteady(params)
    karpata_model2 <- reefSteady(params)
    karpata_model3 <- reefSteady(params)
    
    
    save(karpata_model,   file = "data/karpata_model.rda")
    save(karpata_model2,   file = "data/karpata_model2.rda")
    save(karpata_model3,   file = "data/karpata_model3.rda")
    
    
    
    rm(params)
    
    
    
    
    
   