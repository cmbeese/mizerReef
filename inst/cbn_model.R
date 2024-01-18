# Setting up a Caribbean coral reef mizer model with multiple resources
# Model calibration

#### Setup - loading packages --------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)

#### Parameters ----------------------------------------------------------------

# When new parameters need to be loaded 
setwd("C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/c4_files")

# Load species parameter data
cbn_groups <- read.csv("data/cbn_groups.csv")
cbn_int    <- read.csv("data/cbn_int.csv",  row.names = 1)

## Set model -------------------------------------------------------------------
cbn_model <- newReefParams(group_params = cbn_groups,
                           interaction = cbn_int,
                           method = "noncomplex")

# Finding first steady state ---------------------------------------------------
customFunction("scaleModel", reefScaleModel)
cbn_model <- cbn_model|>
    reef_steady() |> reef_steady() |> reef_steady() 
cbn_model <- calibrateBiomass(cbn_model)

# Errors here again
customFunction("scaleModel", reefScaleModel)
cbn_model <- mizer::matchBiomasses(cbn_model)
cbn_model <- reef_steady(cbn_model) 


# Steady state iteration -------------------------------------------------------
customFunction("scaleModel", reefScaleModel)

cbn_model <- cbn_model |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() |>
    calibrateBiomass() |> matchBiomasses() |> matchGrowth() |> reef_steady() 

# Tune reproduction ------------------------------------------------------------
cbn_model <- setBevertonHolt(cbn_model, erepro = 0.0001)

# and assign it back to the model 
cbn_model <- setBevertonHolt(cbn_model, reproduction_level = rep_level)
cbn_model <- reef_steady(cbn_model)


# Save data --------------------------------------------------------------------
# Params object
save(cbn_model,  file = "cbn_model.rda")
save(cbn_groups, file = "cbn_groups.rda")
save(cbn_int,    file = "cbn_int.rda")

# Copies for package
    save(cbn_model,  
         file = "C:/UsersL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/cbn_model.rda")
    # CSV Files
    save(cbn_groups, 
         file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/cbn_groups.rda")
    save(cbn_int,     
         file = "C:/Users/DELL/OneDrive - Victoria University of Wellington - STAFF/Thesis/345_mizerReef/mizerReef/data/cbn_int.rda")

