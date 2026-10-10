# Create 3-species Caribbean Reef model
# Data object for mizerReef package
# Model calibration

# Load required libraries
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(here)

# Load species and interaction parameters
species_path <- here("inst/data-csv/caribbean_3_species.csv")
interaction_path <- here("inst/data-csv/caribbean_3_interaction.csv")

# Load refuge parameters (tuning profile + Bonaire hole density data)
refuge_path <- here("inst/data-csv/karpata_refuge.csv")
tuning_path <- here("inst/data-csv/tuning_profile.csv")

# Save model parameters as package data objects
save(caribbean_3_species, file = "data/caribbean_3_species.rda")
save(caribbean_3_interaction, file = "data/caribbean_3_interaction.rda")
save(tuning_profile, file = "data/tuning_profile.rda")

# Step 1: Initialize the model with species and interaction parameters
# Note: objects should be initialized with the binned or
#       sigmoidal methods until biomass is calibrated because
#       the competitive method is density-dependent
params <- newReefParams(species_params = caribbean_3_species,
                        interaction = caribbean_3_interaction,
                        method = "binned",
                        method_params = tuning_profile)

# Step 2: moderate density-dependence, then find an initial steady state
rdi <- rep(0.5, nrow(caribbean_3_species))
names(rdi) <- caribbean_3_species$species
params <- setBevertonHolt(params, reproduction_level = rdi)
params <- params |>
  reefSteady() |> reefSteady() |> reefSteady() |>
  reefSteady() |> reefSteady() |> reefSteady()


# Step 3: calibrate the total model biomass and return to steady
params <- calibrateReefBiomass(params)
params <- reefSteady(params)

# Step 4: iterate through matching species-specific growth & biomass,
#         then return to steady
n_iter <- 6 # Note: the ideal number of iterations depends on the model
for (i in 1:n_iter) {
  params <- params |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady()
}

# Step 5: check diets and feeding level; scale plankton abundance if necessary
plotDiet(params)
plotFeedingLevel(params)

# Predators feed on almost entirely plankton; need to decrease plankton
#   abundance kappa with a scaling factor and re-steady the model
scale_down_factor <- 100
params <- scaleDownPlankton(params, scale_down_factor)
for (i in 1:n_iter) {
  params <- params |>
    reefSteady()
}

# Check diets again
plotDiet(params)