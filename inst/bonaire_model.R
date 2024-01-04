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
    bonaire_refuge <- data.frame(L_refuge = 20, prop_protect = 0.4)


# When we can we use saved .rda files
    bonaire_species <- bonaire_species
    bonaire_int     <- bonaire_int
    bonaire_refuge  <- bonaire_refuge
    method <- c("sigmoidal", "noncomplex")
    
# Current steady state attempt:
    # # Remove feeding level for all groups
         bonaire_species$satiation <- rep(FALSE)
    # # Try increasing rho
    #     scale_rho_a <- 2
    #     scale_rho_d <- 1

# Other things tried:
    # # Try widening predation kernels
        # bonaire_species$sigma <- rep(2)
    # # Reduce the proportion of fish that protected
        # bonaire_refuge$prop_protect <- 0.2
    # # Let predators feed equally from all spectra
        # bonaire_int[1,] <- rep(1)

# Base mizer to look at differences with unstructured resources
# params <- newMultispeciesParams(species_params = bonaire_species,
#                                 interaction = bonaire_int,
#                                 min_w_pp = NA,
#                                 w_pp_cutoff = 1.0,
#                                 n = 3/4, p = 3/4)

## Set model -------------------------------------------------------------------
bonaire_model <- newReefParams(species_params = bonaire_species,
                               interaction = bonaire_int,
                               # scale_rho_a = scale_rho_a,
                               # scale_rho_d = scale_rho_d,
                               method = method[1],
                               method_params = bonaire_refuge)

## Project to steady state
bonaire_model <- reef_steady(bonaire_model)
bonaire_model <- calibrateBiomass(bonaire_model)
    # The mizerParams object that goes up to this step is currently saved as 
    # bonaire_model.rda
# The problem comes in this step when I try to match biomasses
bonaire_model <- matchBiomasses(bonaire_model)
bonaire_model <- reef_steady(bonaire_model)

# Scaling up rho_detritus by 2, 4, 5, 10 - model not converging
# Simulation run did not converge after 1.5 years. Value returned by the 
# distance function was: NA
# Warning messages:
#     1: In projectToSteady(params, distance_func = distanceMaxRelRDI, 
#                           t_per = t_per, : inverts are going extinct.
#       2: In tune_algae_detritus(params) : Detrital production is not high 
#                                           enough to support this abundance of
#                                           detritivores. I will increase 
#                                           external input to meet the 
#                                           consumption rate.


# Steady state iteration -------------------------------------------------------
# Doesn't currently work - errors due to inverts going extinct
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

# CSV Files
    save(bonaire_species, file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_species.rda")
    save(bonaire_int,     file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_int.rda")
    save(bonaire_refuge,  file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/bonaire_refuge.rda")

# Build website ----------------------------------------------------------------
pkgdown::build_site()


