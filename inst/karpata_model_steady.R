# Setting up a Caribbean coral reef mizer model with multiple resources
# Model steady state calibration
# last tuned 21st January 2024

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


## Set model -------------------------------------------------------------------
params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

## Project to first steady state -----------------------------------------------
params <- reefSteady(params)

## Calibrate biomasses and growth ----------------------------------------------

# Match observed species group biomasses
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- reefSteady(params)

# Match observed growth rates
params <- matchReefGrowth(params)
params <- reefSteady(params)

# Check for match with age at maturity
age_mat_observed = karpata_species$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Not bad

# Check biomass match - still way off
plotBiomassVsSpecies(params)

# Iterate to refine biomass
params <- params |>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()

plotBiomassVsSpecies(params) # spot on

# Check match with observed age at maturity
age_mat_observed = karpata_species$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Closer than needed

## Now switch to competitive method --------------------------------------------
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

## Check resulting spectra and tune resources ----------------------------------

# Resource looks low - should match sheldon's spectrum
# looks fairly straight not bad but some bumps
plotSpectra(params, total = TRUE, power = 1)
plotSpectra(params, total = TRUE, power = 2)

# plot feeding level to check if resource is too low
plotFeedingLevel(params, species = "inverts")

# Invert feeding level is relatively stable through life, non-linearities
#   are probably due to refuge

# Tune reproduction ------------------------------------------------------------
# We do not have yield or catch data - can't tune size distribution
# First attempt to set very low to see what the minimum values are
params <- setBevertonHolt(params, erepro = 0.0001)
# Now set setting erepro same for all species, as low as possible
params <- setBevertonHolt(params, erepro = 0.04)
# Check reproduction level (value between 0 and 1) - should be higher for
# larger, slow growing species and low for small, fast growing ones
rep <- getReproductionLevel(params)
# These are low for predators, and strangely high for 
# inverts. A reproduction level closer to one means reproduction rate is 
# almost totally independent of the investment into reproduction
# Reproduction should be density independent on reefs

# Check comparison of density dependent & independent reduction
getRDI(params) / getRDD(params)
# Reproduction is density independent for inverts, density dependent for 
# preds but possible too much - RDI is only slightly higher

# Let's increase reproduction level to 0.5 for predators and herbivores
# so that 
rep_level <- c(0.5, rep[2], rep[3])
names(rep_level) <- c("predators","herbivores","inverts")
params <- setBevertonHolt(params,
                          reproduction_level = rep_level)

# Iterate to get back to steady state
params <- params |>
    reefSteady()|>
    reefSteady()|>
    reefSteady()

# Check new reproduction - these look better
rep <- getReproductionLevel(params)
getRDI(params) / getRDD(params)

# Check new spectra
plotSpectra(params, total = TRUE, power = 1)
plotSpectra(params, total = TRUE, power = 2)

# Save!
karpata_model <- reefSteady(params)

# Meh

# Plots ------------------------------------------------------------------------
plotBiomassVsSpecies(karpata_model)
plotRefuge(karpata_model)
plotSpectra(karpata_model, power = 1, total = TRUE)
plotDiet(karpata_model)  
plotGrowthCurves(karpata_model)

# Save in package --------------------------------------------------------------
# Params object
save(karpata_model,   file = "data/karpata_model.rda")

# CSV Files
save(karpata_species, file = "data/karpata_species.rda")
save(karpata_int,     file = "data/karpata_int.rda")
