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
tuning_profile$prop_protect[1:4] <- 0.4

# Attempt 1 --------------------------------------------------------------------
    karpata_10plus$k_vb     <- NULL
    karpata_10plus$l_mat    <- NULL
    # karpata_10plus$alpha <- NULL
    # karpata_10plus$interaction_resource[karpata_10plus$species == "pred_inv"] <- 0.3
    # karpata_10plus$beta[karpata_10plus$species == "pred_inv"] <- 30
    # karpata_10plus$sigma[karpata_10plus$species == "pred_eng"] <- 1.5
    # karpata_10plus$biomass_observed[karpata_10plus$species == "pred_inv"] <- NA
    # karpata_10plus$biomass_cutoff[karpata_10plus$species == "pred_inv"] <- NA
    # karpata_int$pred_inv[karpata_int$pred_inv == 1] <- 0.5
    # karpata_int$pred_inv[karpata_int$pred_inv == 1] <- 0
    
    params <- newReefParams(group_params = karpata_10plus,
                            interaction = karpata_int,
                            method = "binned", w_pp_cutoff = 1,
                            crit_feed = 0.85,
                            method_params = tuning_profile)
    
    # species_params(params)["pred_inv", "pred_kernel_type"] <- "box"
    # species_params(params)["pred_inv", "ppmr_min"] <- 2
    # species_params(params)["pred_inv", "beta"] <- 10
    # species_params(params)["pred_eng", "sigma"] <- 2
    # species_params(params)["pred_grab", "sigma"] <- 2
    # species_params(params)["eels", "sigma"] <- 3
    # species_params(params)["pred_grab", "beta"] <- 80
    
    rdi <- rep(0.8, dim(karpata_int)[1])

    params <- setBevertonHolt(params, reproduction_level = rdi)
    getReproductionLevel(params)
    
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
    # species_params(params)["pred_eng", "sigma"] <- 3
    # params <- matchReefGrowth(params)
    # species_params(params)["pred_grab", "sigma"] <- 3
    # params <- matchReefGrowth(params)
    # species_params(params)["eels", "beta"] <- 100
    # params <- matchReefGrowth(params)
    # species_params(params)["pred_inv", "sigma"] <- 3
    # params <- matchReefGrowth(params)
    
    # betas <- seq(20,400,30)
    # species_params(params)["pred_grab", "beta"] <- betas[1]
    # params <- matchReefGrowth(params)
    # species_params(params)["pred_grab", "beta"] <- betas[2]
    # params <- matchReefGrowth(params)
    # species_params(params)["pred_grab", "beta"] <- betas[3]
    # params <- matchReefGrowth(params)
    # species_params(params)["pred_grab", "beta"] <- betas[4]
    # params <- matchReefGrowth(params)
    # species_params(params)["pred_grab", "beta"] <- betas[5]
    # params <- matchReefGrowth(params)
    # species_params(params)["pred_grab", "beta"] <- betas[6]
    # params <- matchReefGrowth(params)
    # species_params(params)["pred_grab", "beta"] <- betas[7]
    # params <- matchReefGrowth(params)
    # species_params(params)["pred_grab", "beta"] <- betas[8]
    # params <- matchReefGrowth(params)
    
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
    
    # Set to same for all species
    params <- setBevertonHolt(params, erepro = 0.3)
    getReproductionLevel(params)
    getRDD(params)/getRDI(params)
    
    #params <- setBevertonHolt(params, reproduction_level = rdi)
    
    params <- reefSteady(params)
    
    params <- params |>
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
    # plotDiet(karpata_model2)+ scale_x_log10(limits = c(1, 10000))
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
    
    
    
    
    
   