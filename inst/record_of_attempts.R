# Model calibration - record of attempts
# Updated 7th January 2024, new maturity params

# Setup - loading packages -----------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)

# Parameters -------------------------------------------------------------------

# When we can we use saved .rda files
    # as suggested, set L_mat for inverts to same as for preds and herbs
bonaire_species <- bonaire_species
bonaire_int     <- bonaire_int
bonaire_refuge  <- bonaire_refuge
method <- c("sigmoidal", "noncomplex")

# Set model --------------------------------------------------------------------
bonaire_model <- newReefParams(species_params = bonaire_species,
                               interaction = bonaire_int,
                               # scale_rho_a = scale_rho_a,
                               # scale_rho_d = scale_rho_d,
                               # method = method[2])
                               method = method[1],
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
    
    # We need to tell mizer to use `reef_scaleModel()` so that also the
    # consumption parameters are rescaled correctly by `calibrateBiomass()`
    customFunction("scaleModel", reef_scaleModel)
    attempt_1 <- calibrateBiomass(attempt_1)
    attempt_1 <- matchBiomasses(attempt_1)
    attempt_1 <- reef_steady(attempt_1)
    # Convergence is achieved in 12 years but the predators and inverts aren't
    # doing well enough yet, needing an reproductive efficiency greater than 1.

    
    # Simulation run did not converge after 99 years. 
    # Value returned by the distance function was: 689.726402638765
    # Warning message:
    #     In mizer::setBevertonHolt(params, 
    #               reproduction_level = old_reproduction_level) :
    #     The following species require an unrealistic reproductive efficiency 
    #     greater than 1: predators, herbivores, inverts
    
    # Iteration
    attempt_1 <- attempt_1 |>
        calibrateBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Error in steadySingleSpecies(params, species = sel) : 
    #     With these parameter values the herbivores does not have enough food 
    #     to cover its metabolic cost
    
    # getEReproAndGrowth(attempt_1) returns negative values at large sizes - 
    # should it be returning negative values?

## Attempt 2 - change to distanceMaxRelRDI distance function -------------------
    attempt_2 <- reef_steady(bonaire_model, distance_func = distanceMaxRelRDI)
    attempt_2 <- calibrateBiomass(attempt_2)
    attempt_2 <- matchBiomasses(attempt_2)
    # Warning - inverts have unrealistic erepro
    attempt_2 <- reef_steady(attempt_2, distance_func = distanceMaxRelRDI)
    # Error in if (success == TRUE) { : missing value where TRUE/FALSE needed
    attempt_2 <- attempt_2 |>
        calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> 
        reef_steady(distance_func = distanceMaxRelRDI) |>
        calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> 
        reef_steady(distance_func = distanceMaxRelRDI) |>
        calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> 
        reef_steady(distance_func = distanceMaxRelRDI) 
    
    # Error in steadySingleSpecies(params, species = sel) : 
    #       With these parameter values the herbivores does not have enough 
    #       food to cover its metabolic cost

## Attempt 3 - remove feeding level for inverts --------------------------------
    attempt_3 <- bonaire_model
    attempt_3@species_params$satiation <- rep(FALSE)
    attempt_3 <- reef_steady(attempt_3)
    # Converges on fourth or fifth attempt
    attempt_3 <- attempt_3 |>
        reef_steady() |> reef_steady() |> reef_steady() |> reef_steady() 
    # But then errors after calibrating and matching biomasses
    attempt_3 <- attempt_3 |>
        calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
        calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() 
    
    # With these parameter values the herbivores does not have enough food to 
    #   cover its metabolic cost

## Attempt 5 - increase Rho ----------------------------------------------------
    scale_rho_a <- 2
    attempt_5 <- newReefParams(species_params = bonaire_species,
                               interaction = bonaire_int,
                               scale_rho_a = scale_rho_a,
                               # scale_rho_d = scale_rho_d,
                               # method = method[2])
                               method = method[1],
                               method_params = bonaire_refuge)

    attempt_5 <- reef_steady(attempt_5)
    # Converges on fourth or fifth attempt
    attempt_5 <- attempt_5 |>
        reef_steady() |> reef_steady() |> reef_steady()
    # But then errors after calibrating and matching biomasses
    attempt_5 <- attempt_5 |>
        calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
        calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() 
    
    # Error in steadySingleSpecies(params, species = sel) : 
    #       With these parameter values the herbivores does not have enough 
    #       food to cover its metabolic cost
    # In addition: Warning message:
    #     In setBevertonHolt(params) :
    #       For the following species `erepro` has been increased to the smallest 
    #       possible value: erepro[herbivores] = 0.0856; erepro[inverts] = 183000

## Attempt 6 - setting a reasonable invertebrate biomass -----------------------
# is it necessary for me to have data for this?
# AND remove feeding level for them AND increase a
    scale_rho_a <- 4
    attempt_6 <- bonaire_species
    attempt_6$satiation <- rep(FALSE)
    attempt_6$biomass_observed[attempt_6$species == 'inverts'] <- 20
    
    attempt_6 <- newReefParams(species_params = attempt_6,
                                interaction = bonaire_int,
                                scale_rho_a = scale_rho_a,
                                method = method[1],
                                method_params = bonaire_refuge)
    
    # Converges on fourth or fifth attempt
    attempt_6 <- attempt_6 |>
        reef_steady() |> reef_steady() |> reef_steady() |> reef_steady()

    attempt_6 <- attempt_6 |>
        calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
        calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady()
    
    # With scale_rho_a less than 5:
        # Error in steadySingleSpecies(params, species = sel) : 
        #   With these parameter values the herbivores does not have enough 
        #   food to cover its metabolic cost
    
    # With scale_rho_a greater than 5:
        # Error in if (any(growth[idx] == 0)) { : 
        #         missing value where TRUE/FALSE needed
        #  In addition: Warning message:
        #  In setBevertonHolt(params) :
        #       For the following species `erepro` has been increased to the 
        #       smallest possible value: erepro[herbivores] = 493; 
        #       erepro[inverts] = 36300000

# Other things tried: ----------------------------------------------------------

# Getting rid of benthic complexity /predation predation refuge entirely
# nonComplex_bonaire <- newReefParams(species_params = bonaire_species,
#                                       interaction = bonaire_int,
#                                       method = "noncomplex")

# Widening predation kernels
# bonaire_species$sigma <- rep(2)

# Reduce the proportion of fish protected by refuge
# bonaire_refuge$prop_protect <- 0.2

# Scaling up rho_detritus, rho_algae, and both by 2, 4, 5, 10 
# - model not converging
