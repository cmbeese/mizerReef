# Create 3-species Caribbean Reef model
# Data object for mizerReef package
#
# Builds caribbean_3_model from scratch: read the parameter tables, set up the
# model, then project it forward until it settles at steady state.
#
# The model follows the structure of Rogers et al. (2018): herbivores eat only
# algae, invertebrates eat only detritus, and only predators eat plankton.
# The parameter tables already hold the calibrated values, so no calibration
# steps are needed here. They were tuned so that the steady state matches:
#   - biomass above 3.125 g of 107 (predators), 34 (herbivores) and
#     40 (invertebrates) g/m^2
#   - about 10 herbivore juveniles (0.03-1.6 g) per m^2
#   - an adult predator diet of 33% predators, 17% herbivores and
#     50% invertebrates (Rogers et al. 2018)
#   - predators maturing at age 2
#   - algae standing biomass of 350 g/m^2 and a detritus lifetime of 3.85 h
# The calibrated values are: the predators' row of the interaction table;
# w_min, gamma, h, ks, q, erepro and R_max in the species table; herbivores'
# interaction_algae and invertebrates' interaction_detritus, which set the
# algae and detritus scale; and kappa, algae_capacity and external below.

# Load required libraries
library(mizer)
library(mizerReef)
library(here)

# Load species, interaction and refuge parameters
caribbean_3_species <- read.csv(here("inst/data-csv/caribbean_3_species.csv"))
caribbean_3_interaction <- read.csv(
  here("inst/data-csv/caribbean_3_interaction.csv"), row.names = 1)
karpata_refuge <- read.csv(here("inst/data-csv/karpata_refuge.csv"))

# Save model parameters as package data objects
save(caribbean_3_species, file = here("data/caribbean_3_species.rda"))
save(caribbean_3_interaction, file = here("data/caribbean_3_interaction.rda"))

# Step 1: set up the model
params <- newReefParams(
  species_params = caribbean_3_species,
  interaction = caribbean_3_interaction,
  # Refuge: hole densities from Karpata reef in Bonaire; fish compete for
  # holes, and use them from 2.5 cm (0.39 g)
  method = "competitive",
  method_params = karpata_refuge,
  w_settle = 0.025 * 2.5^3,
  # Plankton: individuals up to 10^-1.5 g, the smallest fish size
  w_pp_cutoff = 10^-1.5,
  kappa = 0.1757,
  # Algae has a carrying capacity (g/m^2); detritus does not
  use_UR_cc = TRUE,
  algae_capacity = 902.8,
  detritus_capacity = Inf,
  # External detritus input (g/m^2/year)
  external = 7.757)

# Step 2: project forward from the initial abundances until the model
#         settles at steady state
sim <- project(params, t_max = 300, progress_bar = FALSE)
plotBiomass(sim) # all biomasses should have levelled off
params <- finalParams(sim)

# Step 3: check the steady state against the targets
plotBiomassObservedVsModel(params)
plotDiet(params)
algae_biomass(params)             # target 350 g/m^2
detritus_lifetime(params) * 8760  # target 3.85 hours

# Save the model as a package data object
caribbean_3_model <- params
save(caribbean_3_model, file = here("data/caribbean_3_model.rda"))
