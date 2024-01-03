# Setting up a Caribbean coral reef mizer model with multiple resources
# Model calibration

#### Setup - loading packages and functions ------------------------------------
# Load in relevant packages
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)

#### Parameters ----------------------------------------------------------------

# When new parameters need to be loaded 
    # setwd("C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/c3_files/vignettes")
    # 
    # # Load species parameter data
    # bonaire_species <- read.csv("bonaire_species.csv")
    # bonaire_int     <- read.csv("bonaire_int.csv",  row.names = 1)
    # 
    # # Create some tester refuge scenarios
    # method <- c("sigmoidal", "noncomplex")
    # bonaire_refuge <- data.frame(L_refuge = 20, prop_protect = 0.4)
    # 
# When we can we use saved .rda files
    bonaire_species <- bonaire_species
    bonaire_int     <- bonaire_int
    bonaire_refuge  <- bonaire_refuge
    method <- c("sigmoidal", "noncomplex")
    
# Remove feeding level for now
    bonaire_species$piscivore <- rep(TRUE)

# ## Sorting out unstructured resources
# params <- newMultispeciesParams(species_params = bonaire_species,
#                                 interaction = bonaire_int,
#                                 min_w_pp = NA,
#                                 w_pp_cutoff = 1.0,
#                                 n = 3/4, p = 3/4)

## SET MODEL -------------------------------------------------------------------
bonaire_model <- newReefParams(species_params = bonaire_species,
                               interaction = bonaire_int,
                               scale_rho_a = 1, exp_alg = 0.75,
                               method = method[1],
                               method_params = bonaire_refuge)

## Project to steady state
bonaire_model <- reef_steady(bonaire_model)


bonaire_model <- calibrateBiomass(bonaire_model)
bonaire_model <- matchBiomasses(bonaire_model)

bonaire_model <- bonaire_model |>
calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() 

|>
calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() 

# Plots ------------------------------------------------------------------------
plotBiomassVsSpecies(bonaire_model)
plotRefuge(bonaire_model)
plotSpectra(bonaire_model, power = 1, total = TRUE)
    

# Age at maturity --------------------------------------------------------------
age_mat_observed = bonaire_species$age_mat
age_mat_model = age_mat(bonaire_model)
data.frame(age_mat_model, age_mat_observed)

# Tune feeding level and reproduction using the shiny gadget -------------------
bonaire_model <- tuneParams(bonaire_model)

# Save as rda ------------------------------------------------------------------
save(bonaire_model,   file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_model.rda")
save(bonaire_species, file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_species.rda")
save(bonaire_int,     file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_int.rda")
save(bonaire_refuge,  file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_refuge.rda")




