# Setting up a Caribbean coral reef mizer model with multiple resources
# Model calibration

#### Setup - loading packages --------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)

#### Parameters ----------------------------------------------------------------

# When new parameters need to be loaded 
    setwd("C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/c3_files/vignettes")

    # Load species parameter data
    bonaire_species <- read.csv("bonaire_species.csv")
    bonaire_int     <- read.csv("bonaire_int.csv",  row.names = 1)

    # Create some tester refuge scenarios
    method <- c("sigmoidal", "noncomplex")
    bonaire_refuge <- data.frame(L_refuge = 15, prop_protect = 0.2)

# When we can we use saved .rda files
    # bonaire_species <- bonaire_species
    # bonaire_int     <- bonaire_int
    # bonaire_refuge  <- bonaire_refuge
    # method <- c("sigmoidal", "noncomplex")
    
# Steady state attempts in "inst/record_of_attempts"

# Base mizer to look at differences with unstructured resources
# params <- newMultispeciesParams(species_params = bonaire_species,
#                                 interaction = bonaire_int,
#                                 min_w_pp = NA,
#                                 w_pp_cutoff = 1.0,
#                                 n = 3/4, p = 3/4)

## Set model -------------------------------------------------------------------
bonaire_model <- newReefParams(species_params = bonaire_species,
                               interaction = bonaire_int,
                               method = method[2])
                               # scale_rho_a = scale_rho_a,
                               # scale_rho_d = scale_rho_d,
                               # method = method[1],
                               # method_params = bonaire_refuge)
        
bonaire_model@linecolour["predators"] <-"#8E0408"
bonaire_model@linecolour["inverts"] <- "#D89958"
bonaire_model@linecolour["herbivores"] <- "#578979"
bonaire_model@linecolour["detritus"] <-"burlywood4"
bonaire_model@linecolour["algae"] <- "darkolivegreen4"

## Project to steady state
bonaire_model <- bonaire_model |>
    reef_steady() |> reef_steady() |> reef_steady() |> 
    reef_steady() |> reef_steady() |> reef_steady()

# Steady state iteration -------------------------------------------------------
customFunction("scaleModel", reefScaleModel)
bonaire_model <- bonaire_model |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() 

# Tune reproduction ------------------------------------------------------------
# Do I need to implement fishing to tune resilience?
bonaire_model <- setBevertonHolt(bonaire_model, erepro = 0.0001)

# Save current reproduction level into a vector 
rep_level <- getReproductionLevel(bonaire_model)
# set reproduction level to 0 for all species
rep_level <- c(0,0,0.8)#rep(0)
# and assign it back to the model 
bonaire_model <- setBevertonHolt(bonaire_model, 
                             reproduction_level = rep_level)

bonaire_model <- reef_steady(bonaire_model)

# and plot F curves again
plotYieldVsF(bonaire_model, species = "Predators", 
             F_range = seq(0.1, 0.5, 0.02))

# Cephalopholis cruentata medium resilience species, should have MSY in
# range 0.1 to 0.5 according to mizer course
plotYieldVsF(cel_model, species = "Predators", 
             F_range = seq(0.1, 0.9, 0.02))

# Sparisoma viride medium resilience species, should have MSY in
# range 0.1 to 0.5 according to mizer course
plotYieldVsF(cel_model, species = "Herbivores", 
             F_range = seq(0.1, 0.9, 0.02))

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

    # Params object
    save(bonaire_model,   file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_model.rda")

    
# If I want to save another version of parameters
    # bonaire_model2 <- bonaire_model
    # save(bonaire_model2,   file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_model2.rda")

    # CSV Files
    save(bonaire_species, file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_species.rda")
    save(bonaire_int,     file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_int.rda")
    save(bonaire_refuge,  file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_refuge.rda")

# Build website ----------------------------------------------------------------
pkgdown::build_site()

