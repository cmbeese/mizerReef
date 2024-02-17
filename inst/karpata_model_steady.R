# Setting up a Caribbean coral reef mizer model with multiple resources
# Model steady state calibration
# last tuned 17th February 2024

## Setup - load packages ----------------------------------
library(ggplot2)
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)
library(here)

## Load parameters -----------------------------------
karpata_10plus  <- read.csv(here("inst/karpata_10plus.csv"))
karpata_int     <- read.csv(here("inst/cbn_interaction.csv"),  
                            row.names = 1)
karpata_refuge  <- karpata_refuge
tuning_profile  <- tuning_profile

# Herbivores consume plankton at small sizes and 
#   transition to detritus and algae as they grow
# Invertebrates consume plankton and detritus,
#   with the proportion of detritus increasing with size

## Set model ----------------------------------------
params <- newReefParams(group_params = karpata_10plus,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

## Reduce density dependent of reproduction ----------------
rdi <- rep(0.5, dim(karpata_int)[1])

params <- setBevertonHolt(params, reproduction_level = rdi)
getReproductionLevel(params)

## Project to first steady state -------------------------------
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() 

## Calibrate biomasses and growth ---------------------------------

# Match observed species group biomasses
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- reefSteady(params)

# Match observed growth rates
params <- matchReefGrowth(params)
params <- reefSteady(params)

# Iterate to refine biomass
params <- params |>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()

# Check biomass match
plotBiomassVsSpecies(params) # spot on

# Check match with observed age at maturity
age_mat_observed = karpata_10plus$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Closer than needed

# Check predation mortality, feeding levels, and diets
plotPredMort(params) + facet_wrap(~Species)
plotFeedingLevel(params)
plotDiet(params) + scale_x_log10(limits = c(1, 1e4))
plotSpectra(params, power = 1)

## Now switch to competitive method ----------------------------
params <- newRefuge(params,
                    new_method = "competitive",
                    new_method_params = karpata_refuge)

# Match biomasses again
params <- params |>
    matchBiomasses()|> reefSteady()|> 
    matchBiomasses()|> reefSteady()|>
    matchBiomasses()|> reefSteady()|>
    matchBiomasses()|> reefSteady() 

# Make sure new refuge is in place
plotVulnerable(params)

plotBiomassVsSpecies(params) # spot on

# Check match with observed age at maturity
age_mat_observed = karpata_species$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Still look good

## Check resulting spectra and tune resources------------------------

# Resource looks low - should match sheldon's spectrum
# looks fairly straight not bad but some bumps
plotSpectra(params, total = TRUE, power = 1)
plotSpectra(params, total = TRUE, power = 2)

# plot feeding level to check if resource is too low
plotFeedingLevel(params, species = "inverts")

# Invert feeding level is relatively stable through life, non-linearities
#   are probably due to refuge

# Tune reproduction -----------------------------------------------
# We do not have yield or catch data - can't tune size distribution
# First attempt to set very low to see what the minimum values are
params <- setBevertonHolt(params, erepro = 0.0001)
# Now set setting erepro same for all species, as low as possible
params <- setBevertonHolt(params, erepro = 0.35)
# Project back to steady
params <- reefSteady(params)
# Check reproduction level (value between 0 and 1) - should be higher for
# larger, slow growing species and low for small, fast growing ones
getReproductionLevel(params)
# A reproduction level closer to one means reproduction rate is 
# almost totally independent of the investment into reproduction
# These are near one for all species except farming damsels

# Check comparison of density dependent & independent reduction
getRDI(params) / getRDD(params)
# Reproduction is density independent for almost all species

# Let's increase reproduction level back to 0.5 so there is still some
# density dependence
params <- setBevertonHolt(params, reproduction_level = rdi)

# Iterate to get back to steady state
params <- params |>
    reefSteady()|>
    reefSteady()|>
    reefSteady()

# Check new reproduction - these look better
rep <- getReproductionLevel(params)
getRDI(params) / getRDD(params)

# Check new spectra & plots
plotSpectra(params, total = TRUE, power = 2)
plotPredMort(params) + facet_wrap(~Species)
plotFeedingLevel(params)
plotDiet(params) + scale_x_log10(limits = c(1, 1e4))
plotSpectra(params, power = 1)

# Save!
karpata_model <- reefSteady(params)

# Plots ------------------------------------------
plotBiomassVsSpecies(karpata_model)
plotRefuge(karpata_model)
plotSpectra(karpata_model, power = 1, total = TRUE)
plotDiet(karpata_model)  
plotGrowthCurves(karpata_model)

# Save in package --------------------------------------------
# Params object
save(karpata_model,   file = "data/karpata_model.rda")

# CSV Files
save(karpata_species, file = "data/karpata_species.rda")
save(karpata_int,     file = "data/karpata_int.rda")
