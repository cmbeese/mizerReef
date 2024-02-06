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

## Load parameters --------------------------------------------------

# Load species parameter data
karpata_10plus  <- read.csv(here("inst/karpata_10plus.csv"))
karpata_int     <- read.csv(here("inst/cbn_interaction.csv"),  row.names = 1)
karpata_refuge  <- karpata_refuge
tuning_profile  <- tuning_profile

# Attempt 1 --------------------------------------------------------------------
    karpata_10plus$k_vb     <- NULL
    
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", w_pp_cutoff = 1,
                            crit_feed = 0.85,
                            method_params = tuning_profile,
                            include_ext_mort = FALSE,
                            include_sen_mort = FALSE)
  
    rdi <- rep(0.5, dim(karpata_int)[1])

    params <- setBevertonHolt(params, reproduction_level = rdi)
    getReproductionLevel(params)
    
    ## Project to first steady state
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() 
    
    # Match biomasses
    params <- calibrateReefBiomass(params)
    params <- matchBiomasses(params)
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() 
    
    plotBiomassVsSpecies(params)
    plotSpectra(params, power = 2)
    
    # Now match growth
    params <- matchReefGrowth(params)
    species_params(params)["pred_eng","sigma"] <- 2.5
    species_params(params)["pred_grab","sigma"] <- 3
    species_params(params)["pred_inv","sigma"] <- 3.5
    species_params(params)["eels","sigma"] <- 4
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
    
    plotSpectra(params, power = 2, total = TRUE)
    plotFeedingLevel(params)
    
    # # # Scale down background??
    params <- scaleReefBackground(params, 2)
    # 
    params <- params |>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady()|>
        calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
        reefSteady() |> reefSteady()

    # Check for match with age at maturity
    age_mat_observed = karpata_10plus$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)
    # Check biomass match
    plotBiomassVsSpecies(params)

    rep_level <- getReproductionLevel(params)
    # then we replace our species' reproduction 
    #level with a new value 
    rep_level["farm_damsel"] <- 0.01
    # and assign it back to the model 
    params <- setBevertonHolt(params,
                              reproduction_level = rep_level)
    
    
    params <- setBevertonHolt(params, erepro = 0.0001)
    
    # Set to same for all species
    rdi <- rep(0.5, dim(karpata_int)[1])
    params <- setBevertonHolt(params, erepro = 0.08)
    getReproductionLevel(params)
    
    params <- setBevertonHolt(params, reproduction_level = 0.1)
    
    params <- setBevertonHolt(params, reproduction_level = rdi)
    params <- reefSteady(params)
    params <- params |>
        reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady() |> reefSteady() |> 
        reefSteady() |> reefSteady()
    getReproductionLevel(params)
    
    
    plotSpectra(params, power = 2) + 
        geom_vline(xintercept = 1600) +
        geom_vline(xintercept = 100) +
        geom_vline(xintercept = 50)
    
    plotPredMort(params) + facet_wrap(~Species)
    plotFeedingLevel(params)
    plotDiet(params) + scale_x_log10(limits = c(0.1, 1e4))
    plotDiet(params) + scale_x_log10(limits = c(1, 1e4))
    plotDiet(karpata_model)+ scale_x_log10(limits = c(1, 1e4))
    params <- tuneParams(params)
    
    plotDiet(karpata_model2)+ scale_x_log10(limits = c(1, 1e4),
                                            breaks = c(1,10,100,
                                                       1e3,1e4),
                                            labels = c(1,10,100,
                                                       "1e3","1e4"))
    plotDiet(karpata_model3)+ scale_x_log10(limits = c(1, 1e4),
                                            breaks = c(1,10,100,
                                                       1e3,1e4),
                                            labels = c(1,10,100,
                                                       "1e3","1e4"))
    
## Now switch to competitive method ----------------------------
    params <- newRefuge(params,
                        new_method = "competitive",
                        new_method_params = karpata_refuge)
    
    # Match biomasses again
    params <- params |>
        matchBiomasses()|> reefSteady()|> 
        matchBiomasses()|> reefSteady()|>
        matchBiomasses()|> reefSteady()|>
        matchBiomasses()|> reefSteady()|> 
        reefSteady()|> reefSteady()|> reefSteady()|> 
        reefSteady()|> reefSteady()
    
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
    
    
    
    
    
   