library(mizer)
library(mizerReef)

params <- newReefParams(group_params = karpata_species,
                        interaction = karpata_int,
                        method = "competitive",
                        method_params = karpata_refuge)

params@other_params$degrade = FALSE

params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() 


params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- reefSteady(params)
