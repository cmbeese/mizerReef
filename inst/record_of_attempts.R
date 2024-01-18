# Model calibration - record of attempts
# Updated 18th January 2024, new invertebrate params

# Setup - loading packages -----------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)

# Parameters -------------------------------------------------------------------
bonaire_species <- bonaire_species
bonaire_int     <- bonaire_int
bonaire_refuge  <- bonaire_refuge

# Set model --------------------------------------------------------------------
bonaire_model <- newReefParams(group_params = bonaire_species,
                               interaction = bonaire_int,
                               method = "competitive",
                               method_params = bonaire_refuge)
# We need to specify linecolours for detritus and algae, for them to be included
# in diet plots
bonaire_model@linecolour["detritus"] <-"brown"
bonaire_model@linecolour["algae"] <- "darkseagreen"

# Project to steady state ------------------------------------------------------
# Changed default distance function back to distanceSSLogN, but trying both to
# see if either work

## Attempt 1 - no changes-------------------------------------------------------
    attempt_1 <- reef_steady(bonaire_model)
    
    # Converges on fourth or fifth attempt
    attempt_1 <- attempt_1 |>
        reef_steady() |> reef_steady() |> reef_steady() |> 
        reef_steady() |> reef_steady()
    
    # Iteration
    attempt_1 <- attempt_1 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Error in steadySingleSpecies(params, species = sel) : 
    #     With these parameter values the predators does not have enough food 
    #     to cover its metabolic cost
    
## Attempt 2 - widening the predation kernel -----------------------------------
    attempt_2 <- bonaire_model
    attempt_2 <- set_species_param_default(attempt_2, "sigma", 2)  
    
    attempt_2 <- reef_steady(attempt_2)
    
    # Converges on fourth or fifth attempt
    attempt_2 <- attempt_2 |>
        reef_steady() |> reef_steady() |> reef_steady() |> 
        reef_steady() |> reef_steady()
    
    # Iteration
    attempt_2 <- attempt_2 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Error in steadySingleSpecies(params, species = sel) : 
    #     With these parameter values the predators does not have enough food 
    #     to cover its metabolic cost
    
    # In addition: Warning message:
    #     In setBevertonHolt(params) :
    #     For the following species `erepro` has been increased to the 
    #       smallest possible value: erepro[herbivores] = 4.48e-05; 
    #       erepro[inverts] = 0.0636
    
## Attempt 3 - widening the predation kernel, setting all interactions to 1 ----
    attempt_3 <- bonaire_model
    attempt_3 <- set_species_param_default(attempt_3, "sigma", 2)  
    attempt_3@interaction["predators", ] <- rep(1)
    attempt_3 <- reef_steady(attempt_3)
    
    # Converges immediately!
    attempt_3 <- calibrateReefBiomass(attempt_3)
    attempt_3 <- matchBiomasses(attempt_3)
    
    # In setBevertonHolt(params) :
    #     For the following species `erepro` has been increased to the smallest possible value: erepro[predators] = 0.0528; erepro[herbivores] = 0.000216; erepro[inverts] = 0.0263
    
    attempt_3 <- matchGrowth(attempt_3)
    
    # This is the step that it doesn't like
    
    # Error in steadySingleSpecies(params, species = sel) : 
    #     With these parameter values the predators does not have enough food 
    #     to cover its metabolic cost
    
## Attempt 4 - widening the predation kernel, setting all interactions to 1,
    # scale up background resource ----
    attempt_4 <- bonaire_model
    attempt_4 <- set_species_param_default(attempt_4, "sigma", 2)  
    attempt_4@interaction["predators", ] <- rep(1)
    attempt_4 <- reef_steady(attempt_4)
    attempt_4 <- scaleDownBackground(attempt_4, factor = 1/2)
    
    # Converges immediately!
    attempt_4 <- calibrateReefBiomass(attempt_4)
    attempt_4 <- matchBiomasses(attempt_4)
    attempt_4 <- reef_steady(attempt_4)
    # Detritus keeps needing to be scaled up?
    attempt_4 <- matchGrowth(attempt_4)
    
    # Iteration
    attempt_4 <- attempt_4 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Error in steadySingleSpecies(params, species = sel) : 
    #     With these parameter values the predators does not have enough food 
    #     to cover its metabolic cost
