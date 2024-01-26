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
karpata_int     <- read.csv(here("inst/cbn_interaction.csv"),  row.names = 1)
karpata_refuge  <- karpata_refuge
tuning_profile  <- tuning_profile

# With these parameters, herbivores consume plankton at small sizes and 
#   transition fully to algae by maturity 
# With these parameters, invertebrates consume plankton and detritus,
#   with the proportion of detritus increasing with size

# Attempt 1 --------------------------------------------------------------------
## Set model
params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

## Project to first steady state
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
# Converges on 7th attempt

## Calibrate biomasses and growth
# Match observed species group biomasses
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

plotSpectra(params)

params <- setBevertonHolt(params, erepro = 0.0001)   

# Match observed growth rates
params <- matchReefGrowth(params)
params <- reefSteady(params)

# Error in steadySingleSpecies(params, species = sel) : 
#     With these parameter values the pred_grab does not have enough
#     food to cover its metabolic cost

# Check for match with age at maturity
age_mat_observed = karpata_species$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Terrible

rm(params)

# Attempt 2 --------------------------------------------------------------------
# Widening the predation kernel for grabbing predators
karpata_species$sigma[karpata_species$species == "pred_grab"] <- 2
    
params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

## Project to first steady state
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
# Converges on 7th attempt

## Calibrate biomasses and growth
# Match observed species group biomasses
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

plotSpectra(params)

# Match observed growth rates
params <- matchReefGrowth(params)
params <- reefSteady(params)

# Error in steadySingleSpecies(params, species = sel) : 
#     With these parameter values the pred_grab does not have enough food to 
#     cover its metabolic cost

rm(params)

# Attempt 3 --------------------------------------------------------------------
# Setting resource interaction to 1 for detritivores and herbivores
# also widened predation kernel for grabbing predators
karpata_species$interaction_algae[karpata_species$interaction_algae == 0.5] <- 1
karpata_species$interaction_algae[karpata_species$interaction_detritus == 0.5] <- 1
karpata_species$sigma[karpata_species$species == "pred_grab"] <- 2

params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

## Project to first steady state
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
# Converges on 7th attempt

## Calibrate biomasses and growth
# Match observed species group biomasses
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

# Match observed growth rates
params <- matchReefGrowth(params)
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() 

# Check for match with age at maturity
age_mat_observed = karpata_species$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Very close

# Check biomass match - still way off
plotBiomassVsSpecies(params)

# Iterate to refine biomass
params <- params |>
    calibrateReefBiomass() |> matchBiomasses() |> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    matchReefGrowth()|> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses() |> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    matchReefGrowth()|> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses() |> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady()

# Error in steadySingleSpecies(params, species = sel) : 
#     With these parameter values the pred_grab does not have enough food to 
#     cover its metabolic cost

plotBiomassVsSpecies(params) # still way off :()

rm(params)

# Attempt 4 --------------------------------------------------------------------
# Setting resource interaction to 1 for detritivores and herbivores
# also widened predation kernel for grabbing predators
# Same as before but matching biomasses species by species
karpata_species$interaction_algae[karpata_species$interaction_algae == 0.5] <- 1
karpata_species$interaction_algae[karpata_species$interaction_detritus == 0.5] <- 1
karpata_species$sigma[karpata_species$species == "pred_grab"] <- 2

params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

## Project to first steady state
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
# Converges on 7th attempt

## Calibrate biomasses and growth
# Match observed species group biomasses
params <- calibrateReefBiomass(params)
# Match biomass species by species
# Predators first
params <- matchBiomasses(params, species = "pred_eng")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_grab")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_plank")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_inv")
params <- reefSteady(params)

# Now herbivores
params <- matchBiomasses(params, species = "herbs")
params <- reefSteady(params)
# Stops converging here - does order matter?

params <- matchBiomasses(params, species = "farm_damsel")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "parrotfish")
params <- reefSteady(params)

params <- matchBiomasses(params)
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |> 
    reefSteady() |> reefSteady() |> reefSteady()

# Match observed growth rates
params <- matchReefGrowth(params)
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() 

# Check for match with age at maturity
age_mat_observed = karpata_species$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Very close

# Check biomass match - still way off
plotBiomassVsSpecies(params)

# Iterate to refine biomass
params <- params |>
    calibrateReefBiomass() |> matchBiomasses() |> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    matchReefGrowth()|> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses() |> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    matchReefGrowth()|> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses() |> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses() |> 
    reefSteady() |> reefSteady() |> reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses() |> 
    reefSteady() |> reefSteady() |> reefSteady()

# Error in mizer::setBevertonHolt(params, reproduction_level = old_reproduction_level) : 
#     Some species have no reproduction.
plotBiomassVsSpecies(params) # still way off :()

rm(params)


# Attempt 4 --------------------------------------------------------------------
# Setting resource interaction to 1 for detritivores and herbivores
# also widened predation kernel for grabbing predators
# Matching biomasses species by species, herbivores first and ensuring
# convergence between each species
karpata_species$interaction_algae[karpata_species$interaction_algae == 0.5] <- 1
karpata_species$interaction_algae[karpata_species$interaction_detritus == 0.5] <- 1
karpata_species$sigma[karpata_species$species == "pred_grab"] <- 2

params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

## Project to first steady state
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
# Converges on 7th attempt

## Calibrate biomasses and growth
# Match observed species group biomasses
params <- calibrateReefBiomass(params)
# Match biomass species by species
# herbivores
params <- matchBiomasses(params, species = "herbs")
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

params <- matchBiomasses(params, species = "farm_damsel")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "parrotfish")
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

# Predators
params <- matchBiomasses(params, species = "pred_inv")
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

params <- matchBiomasses(params, species = "pred_plank")
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

params <- matchBiomasses(params, species = "pred_eng")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_grab")
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

# Match all 
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- reefSteady(params)

# Check biomass match - spot on 
plotBiomassVsSpecies(params)

# Match observed growth rates
params <- matchReefGrowth(params)

# Error in steadySingleSpecies(params, species = sel) : 
#     With these parameter values the pred_grab does not have enough food to 
#     cover its metabolic cost

# Can't do iteration
rm(params)

# Attempt 5 --------------------------------------------------------------------
# Increase refuge proportion in the tuning profile
# Setting resource interaction to 1 for detritivores and herbivores
# also widened predation kernel for grabbing predators
# Matching biomasses species by species, herbivores first and ensuring
# convergence between each species
tuning_profile$prop_protect <- 2*tuning_profile$prop_protect

karpata_species$interaction_algae[karpata_species$interaction_algae == 0.5] <- 1
karpata_species$interaction_algae[karpata_species$interaction_detritus == 0.5] <- 1
karpata_species$sigma[karpata_species$species == "pred_grab"] <- 2

params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

## Project to first steady state
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
# Converges on 7th attempt

## Calibrate biomasses and growth
# Match observed species group biomasses
params <- calibrateReefBiomass(params)
# Match biomass species by species
# Easy species first
params <- matchBiomasses(params, species = "pred_inv")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_plank")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_eng")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_grab")
params <- reefSteady(params)

# herbivores
params <- matchBiomasses(params, species = "farm_damsel")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "parrotfish")
params <- reefSteady(params)

plotBiomassVsSpecies(params)

# Herbivores are the problem species
params <- matchBiomasses(params, species = "herbs")
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

# Match all 
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()

# Check biomass match - spot on 
plotBiomassVsSpecies(params)

# Match observed growth rates
params <- matchReefGrowth(params)

# Error in steadySingleSpecies(params, species = sel) : 
#     With these parameter values the pred_grab does not have enough food to 
#     cover its metabolic cost

# Can't do iteration
rm(params)

# Attempt 6 --------------------------------------------------------------------
# Increase refuge proportion in the tuning profile
# Let mizer set the default biomass for herbivores
# Setting resource interaction to 1 for detritivores and herbivores
# also widened predation kernel for grabbing predators
# Matching biomasses species by species, herbivores first and ensuring
# convergence between each species
tuning_profile$prop_protect <- 2*tuning_profile$prop_protect
karpata_species$biomass_observed[karpata_species$species == "herbs"] <- NA
karpata_species$biomass_cutoff[karpata_species$species == "herbs"] <- NA

karpata_species$interaction_algae[karpata_species$interaction_algae == 0.5] <- 1
karpata_species$interaction_algae[karpata_species$interaction_detritus == 0.5] <- 1
karpata_species$sigma[karpata_species$species == "pred_grab"] <- 2

params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

## Project to first steady state
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady()
# Converges on 7th attempt

## Calibrate biomasses and growth
# Match observed species group biomasses
params <- calibrateReefBiomass(params)
# Match biomass species by species
# Easy species first
params <- matchBiomasses(params, species = "pred_inv")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_plank")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_eng")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_grab")
params <- reefSteady(params)

# herbivores
params <- matchBiomasses(params, species = "farm_damsel")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "parrotfish")
params <- reefSteady(params)

plotBiomassVsSpecies(params)

# Match all 
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- reefSteady(params)

# Check biomass match - spot on 
plotBiomassVsSpecies(params)

# Bumpy
plotSpectra(params)

# Match observed growth rates
params <- matchReefGrowth(params)

# Error in steadySingleSpecies(params, species = sel) : 
#     With these parameter values the pred_grab does not have enough food to 
#     cover its metabolic cost

# Can't do iteration
rm(params)

# Attempt 7 --------------------------------------------------------------------
# Reduced w_pp_cutoff size - resource shouldnt overlap so much with fish spectra
# Set to same as w_settle
# Same changes as attempt 6
params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile,
                        w_pp_cutoff = 0.1)

## Project to first steady state
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() 

## Calibrate biomasses and growth
# Match observed species group biomasses
params <- calibrateReefBiomass(params)
# Match biomass species by species
# Easy species first
params <- matchBiomasses(params, species = "pred_inv")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_plank")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_eng")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "pred_grab")
params <- reefSteady(params)

# herbivores
params <- matchBiomasses(params, species = "farm_damsel")
params <- reefSteady(params)

params <- matchBiomasses(params, species = "parrotfish")
params <- reefSteady(params)

plotBiomassVsSpecies(params)

# Match all 
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- reefSteady(params)

# Check biomass match - spot on 
plotBiomassVsSpecies(params)

# Bumpy
plotSpectra(params)

# Match observed growth rates
params <- matchReefGrowth(params)

# Check biomass match - spot on 
plotBiomassVsSpecies(params)

# Iterate to refine biomass
params <- params |>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()

# Error in steadySingleSpecies(params, species = sel) : 
#     With these parameter values the pred_grab does not have enough food to 
#     cover its metabolic cost
