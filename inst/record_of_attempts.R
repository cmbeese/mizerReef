# Model calibration - record of attempts

# Setup - loading packages -----------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)

# Parameters -------------------------------------------------------------------

# When we can we use saved .rda files
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

bonaire_model@linecolour["predators"] <-"#8E0408"
bonaire_model@linecolour["inverts"] <- "#D89958"
bonaire_model@linecolour["herbivores"] <- "#578979"

# Project to steady state ------------------------------------------------------
# Changed default distance function back to distanceSSLogN, but trying both to
# see if either work

## Attempt 1 - no changes-------------------------------------------------------
attempt_1 <- reef_steady(bonaire_model)

# Converges on fourth or fifth attempt
attempt_1 <- attempt_1 |>
    reef_steady() |> reef_steady() |> reef_steady() |> 
    reef_steady() |> reef_steady() |> reef_steady()

# But then errors after calibrating and matching biomasses
attempt_1 <- calibrateBiomass(attempt_1)
attempt_1 <- matchBiomasses(attempt_1)
attempt_1 <- reef_steady(attempt_1)

# Simulation run did not converge after 99 years. 
# Value returned by the distance function was: 630.781128743289
# Error in mizer::setBevertonHolt(params, 
#                               reproduction_level = old_reproduction_level): 
# Some species have no reproduction.

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
#       With these parameter values the herbivores does not have enough food to 
#       cover its metabolic cost

## Attempt 3 - remove feeding level for inverts --------------------------------
attempt_3 <- bonaire_model
attempt_3@species_params$satiation <- rep(FALSE)
attempt_3 <- reef_steady(attempt_3)
# Converges on fourth or fifth attempt
attempt_3 <- attempt_3 |>
    reef_steady() |> reef_steady() |> reef_steady() |> 
    reef_steady() |> reef_steady() |> reef_steady()
# But then errors after calibrating and matching biomasses
attempt_3 <- attempt_3 |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() 

# Error in setBevertonHolt(params) : Some species have no reproduction.

## Attempt 4 - same as attempt 3 with distanceMaxRelRDI distance function ------
attempt_4 <- bonaire_model
attempt_4@species_params$satiation <- rep(FALSE)
attempt_4 <- reef_steady(attempt_4, distance_func = distanceMaxRelRDI)
attempt_4 <- attempt_4 |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> 
    reef_steady(distance_func = distanceMaxRelRDI) |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> 
    reef_steady(distance_func = distanceMaxRelRDI) |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> 
    reef_steady(distance_func = distanceMaxRelRDI) 

# Error in steadySingleSpecies(params, species = sel) :
#       With these parameter values the herbivores does not have enough food to 
#       cover its metabolic cost
# In addition: Warning message:
#     In setBevertonHolt(params) :
#       For the following species `erepro` has been increased to the smallest 
#       possible value: erepro[herbivores] = 0.00058; erepro[inverts] = 3.22e+08

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
    reef_steady() |> reef_steady() |> reef_steady() |> 
    reef_steady() |> reef_steady() |> reef_steady()
# But then errors after calibrating and matching biomasses
attempt_5 <- attempt_5 |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() 

# Unrealistic reproduction for invertebratesssssssssss
# Error in setBevertonHolt(params) : Some species have no reproduction.

## Attempt 6 - setting a reasonable invertebrate biomass -----------------------
# is it necessary for me to have data for this?
bonaire_species$biomass_observed[bonaire_species$species == 'inverts'] <- 20
attempt_6 <- newReefParams(species_params = bonaire_species,
                           interaction = bonaire_int,
                           scale_rho_a = scale_rho_a,
                           method = method[1],
                           method_params = bonaire_refuge)

attempt_6 <- reef_steady(attempt_6)
attempt_6 <- attempt_6 |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady()

# Error in setBevertonHolt(params) : Some species have no reproduction.

# Other things tried: ----------------------------------------------------------

# Getting rid of betnhic complexity /predation predation refuge entirely
# nonComplex_bonaire <- newReefParams(species_params = bonaire_species,
#                                       interaction = bonaire_int,
#                                       method = "noncomplex")

# Widening predation kernels
# bonaire_species$sigma <- rep(2)

# Reduce the proportion of fish protected by refuge
# bonaire_refuge$prop_protect <- 0.2

# Scaling up rho_detritus, rho_algae, and both by 2, 4, 5, 10 
# - model not converging