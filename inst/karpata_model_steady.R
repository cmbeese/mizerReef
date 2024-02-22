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
                        w_pp_cutoff = 1,
                        method = "binned",
                        method_params = tuning_profile)

## Reduce density dependent of reproduction ----------------
rdi <- rep(0.5, dim(karpata_int)[1])

params <- setBevertonHolt(params, reproduction_level = rdi)
getReproductionLevel(params)

## Project to first steady state -------------------------------
params <- reefSteady(params)

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
    reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()

# Check biomass match
plotBiomassVsSpecies(params) # spot on
plotTotalAbundance(params)
plotTotalBiomass(params)

# Check match with observed age at maturity
age_mat_observed = karpata_10plus$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Closer than needed

# Check predation mortality, feeding levels, and diets
plotPredMort(params) + facet_wrap(~Species)
plotDiet(params) + scale_x_log10(limits = c(0.1, 1e4),
                                 breaks = c(1, 10, 100, 1000))

## Check resulting spectra and tune resources-----------------------------------

# Resource looks low - should match sheldon's spectrum
# looks fairly straight not bad but some bumps
plotSpectra(params, total = TRUE, power = 1)
plotSpectra(params, total = TRUE, power = 2)

# plot feeding level to check if resource is too low
plotFeedingLevel(params)

params <- scaleReefBackground(params, factor = 1.5)
# Invert feeding level is relatively stable through life, non-linearities
#   are probably due to refuge

# Iterate to refine biomass
params <- params |>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()

# Resource looks low - should match sheldon's spectrum
# looks fairly straight not bad but some bumps
plotSpectra(params, power = 1)
plotSpectra(params, total = TRUE, power = 2)

# plot feeding level to check if resource is too low
plotFeedingLevel(params)


## Now switch to competitive method ----------------------------
params <- newRefuge(params,
                    new_method = "competitive",
                    new_method_params = karpata_refuge)

# Iterate to refine biomass
params <- params |>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|> 
    reefSteady()

# Make sure new refuge is in place
plotVulnerable(params)
plotRefuge(params)
plotBiomassVsSpecies(params) # spot on

# Check match with observed age at maturity
age_mat_observed = karpata_10plus$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Still look good

# Tune reproduction -----------------------------------------------
# We do not have yield or catch data - can't tune size distribution
# First attempt to set very low to see what the minimum values are
params <- setBevertonHolt(params, erepro = 0.0001)
# Now set setting erepro same for all species, as low as possible
params <- setBevertonHolt(params, erepro = 0.0006)
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
test4 <- reefSteady(params)
test_sp4 <- karpata_10plus
test_i4 <- karpata_int

# Save in package --------------------------------------------
save(test4,    file = "data/test4.rda")
save(test_sp4, file = "data/test_sp4.rda")
save(test_i4,  file = "data/test_i4.rda")



# Params object
save(karpata_model3,   file = "data/karpata_model3.rda")

# CSV Files
save(karpata_species, file = "data/karpata_species.rda")
save(karpata_int,     file = "data/karpata_int.rda")
