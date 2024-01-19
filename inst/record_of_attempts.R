# Model calibration - record of attempts
# Updated 19th January

# Setup - loading packages -----------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)
# Parameters -------------------------------------------------------------------
bonaire_species <- bonaire_species
bonaire_int     <- bonaire_int
bonaire_refuge  <- bonaire_refuge
# Set model --------------------------------------------------------------------
bonaire_model <- newReefParams(group_params = bonaire_species,
                               interaction = bonaire_int,
                               method = "competitive",
                               method_params = bonaire_refuge)
# Changes to code since 18th Jan -----------------------------------------------
    # reefScaleModel - added `params@other_params$algae$growth <- 
    # params@other_params$algae$growth * factor` to scale algae growth with
    # calibrate biomass - not sure if this is right/ necessary but I noticed
    # it was missing before, algae was getting very high

## Attempt 1 -------------------------------------------------------------------
# no changes
    attempt_1 <- reef_steady(bonaire_model)
    
    # Converges on fourth or fifth attempt
    attempt_1 <- attempt_1 |>
        reef_steady() |> reef_steady() |> reef_steady() |> 
        reef_steady() |> reef_steady()
    
    # Calibrate
    attempt_1 <- calibrateReefBiomass(attempt_1)
    attempt_1 <- matchBiomasses(attempt_1)
    attempt_1 <- reef_steady(attempt_1)
    
    # Growth
    age_mat_model = age_mat(attempt_1)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Match growth rates
    attempt_1 <- matchGrowth(attempt_1)
    
    # Growth - still way off
    age_mat_model = age_mat(attempt_1)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Iteration
    attempt_1 <- attempt_1 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    # Predators and inverts need reproductive efficiency greater than 1
    
    # Look at spectra - look crazy!
    plotSpectra(attempt_1, power = 2)
    
    # Increase resource
    attempt_1 <- scaleDownBackground(attempt_1, factor = 1/2)
    
    # Find steady again
    attempt_1 <- attempt_1 |> matchGrowth() |> reef_steady()
    
    # Maturation age for herbivores still very low
    age_mat_model = age_mat(attempt_1)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Increase resource again
    attempt_1 <- scaleDownBackground(attempt_1, factor = 1/2)
    
    # Find steady again
    attempt_1 <- attempt_1 |> matchGrowth() |> reef_steady()
    
    # Reproduction
    attempt_1 <- setBevertonHolt(attempt_1, erepro = 0.0001)
    # For the following species `erepro` has been increased to the smallest 
    # possible value: erepro[predators] = 349; erepro[inverts] = 3.01
    
    plotSpectra(attempt_1, power = 2)
    # Herbivores seem to be going crazy?
    
    rm(attempt_1)
    
## Attempt 2 -------------------------------------------------------------------
    # widening the predation kernel 
    attempt_2 <- bonaire_model
    attempt_2 <- set_species_param_default(attempt_2, "sigma", 2)  
    attempt_2 <- reef_steady(attempt_2)
    
    # Converges on fourth or fifth attempt
    attempt_2 <- attempt_2 |>
        reef_steady() |> reef_steady() |> reef_steady() |> 
        reef_steady() |> reef_steady()
    
    # Calibrate
    attempt_2 <- calibrateReefBiomass(attempt_2)
    attempt_2 <- matchBiomasses(attempt_2)
    attempt_2 <- reef_steady(attempt_2)
    
    # Growth
    age_mat_model = age_mat(attempt_2)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Match growth rates
    attempt_2 <- matchGrowth(attempt_2)
    attempt_2 <- reef_steady(attempt_2)
    
    # Growth - age mat for preds goes way up here?
    age_mat_model = age_mat(attempt_2)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Iteration
    attempt_2 <- attempt_2 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Fails on third line 
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 64.8547881807725
    
    # If I run it again: With these parameter values the predators does not 
    # have enough food to cover its metabolic cost
    
    rm(attempt_2)

## Attempt 3 -------------------------------------------------------------------
#  widening the predation kernel, setting all interactions to 1
    attempt_3 <- bonaire_model
    attempt_3 <- set_species_param_default(attempt_3, "sigma", 2)  
    attempt_3@interaction["predators", ] <- rep(1)
    attempt_3 <- reef_steady(attempt_3)
    
    # Converges immediately!
    attempt_3 <- calibrateReefBiomass(attempt_3)
    attempt_3 <- matchBiomasses(attempt_3)
    attempt_3 <- reef_steady(attempt_3)
    
    # Growth - age mat low for herbs
    age_mat_model = age_mat(attempt_3)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Match growth rates
    attempt_3 <- matchGrowth(attempt_3)
    attempt_3 <- reef_steady(attempt_3)
    
    # Growth - age mat still low for herbs
    age_mat_model = age_mat(attempt_3)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Iteration
    attempt_3 <- attempt_3 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 239.302546379538
    # Error in if (any(growth[idx] == 0)) { : 
    #         missing value where TRUE/FALSE needed
    
    # Distance function returns same value every time
    
    rm(attempt_3)

## Attempt 4 -------------------------------------------------------------------
    # widening the predation kernel, setting all interactions to 1 and 
    # setting gamma to value from Alice's model - 6.4 for preds, 0.2 inv and herb
    attempt_4 <- bonaire_model
    attempt_4 <- set_species_param_default(attempt_4, "sigma", 2)
    attempt_4@species_params$gamma <- c(6.4,0.2,0.2)
    attempt_4@interaction["predators", ] <- rep(1)
    attempt_4 <- reef_steady(attempt_4)
    
    # Converges immediately!
    attempt_4 <- calibrateReefBiomass(attempt_4)
    attempt_4 <- matchBiomasses(attempt_4)
    attempt_4 <- reef_steady(attempt_4)
    
    # Detritus keeps needing to be scaled up?
    attempt_4 <- matchGrowth(attempt_4)
    attempt_4 <- reef_steady(attempt_4)
    
    # Iteration
    attempt_4 <- attempt_4 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 239.302546379538
    # Error in if (any(growth[idx] == 0)) { : 
    #         missing value where TRUE/FALSE needed
    
    # Same exact distance function value as attempt 3 - that's weird right?
    
    rm(attempt_4)
    
## Attempt 5 -------------------------------------------------------------------
    # widening the predation kernel, setting all interactions to 1, and
    # increasing gamma to higher values from range in Nash 2015 
    attempt_5 <- bonaire_model
    attempt_5 <- set_species_param_default(attempt_5, "sigma", 2)  
    attempt_5@interaction["predators", ] <- rep(1)
    attempt_5@species_params$gamma <- c(100,5,0.2)
    attempt_5 <- reef_steady(attempt_5)
    
    # Converges immediately!
    attempt_5 <- calibrateReefBiomass(attempt_5)
    attempt_5 <- matchBiomasses(attempt_5)
    attempt_5 <- reef_steady(attempt_5)
    
    # Detritus keeps needing to be scaled up?
    attempt_5 <- matchGrowth(attempt_5)
    attempt_5 <- reef_steady(attempt_5)
    
    # Iteration
    attempt_5 <- attempt_5 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 239.302546379538
    # Error in if (any(growth[idx] == 0)) { : 
    #         missing value where TRUE/FALSE needed
    
    # Again same distance function return
    
    rm(attempt_5)
    
## Attempt 6 -------------------------------------------------------------------
# widening the predation kernel, setting all interactions to 1, and
# increasing gamma to higher values from range in Nash 2015, let invertebrates 
# get bigger so predators can feed on them as well
    attempt_6 <- bonaire_model
    attempt_6 <- set_species_param_default(attempt_6, "sigma", 2)  
    attempt_6@interaction["predators", ] <- rep(1)
    attempt_6@species_params$l_max <- c(50,50,30)
    attempt_6@species_params$gamma <- c(100,NA,NA)
    attempt_6 <- reef_steady(attempt_6)
    
    # Converges immediately!
    attempt_6 <- calibrateReefBiomass(attempt_6)
    attempt_6 <- matchBiomasses(attempt_6)
    attempt_6 <- reef_steady(attempt_6)
    
    # Growth
    attempt_6 <- matchGrowth(attempt_6)
    attempt_6 <- reef_steady(attempt_6)
    
    
    attempt_6 <- attempt_6 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 239.302546379538
    # Error in if (any(growth[idx] == 0)) { : 
    #         missing value where TRUE/FALSE needed
    
    # SAME DISTANCE AGAIN
    
    rm(attempt_6)
    
## Attempt 7 -------------------------------------------------------------------
    # widening the predation kernel, setting all interactions to 1, and
    # increasing gamma to higher values from range in Nash 2015, let invertebrates 
    # get bigger so predators can feed on them as well, scale up resource
    attempt_7 <- bonaire_model
    attempt_7 <- set_species_param_default(attempt_7, "sigma", 2)  
    attempt_7@interaction["predators", ] <- rep(1)
    attempt_7@species_params$l_max <- c(50,50,30)
    attempt_7@species_params$gamma <- c(100,NA,NA)
    attempt_7 <- scaleDownBackground(attempt_7, factor = 1/2)
    attempt_7 <- reef_steady(attempt_7)
    attempt_7 <- reef_steady(attempt_7)
    attempt_7 <- reef_steady(attempt_7)
    # Doesn't converge - and distance function is the same every time?
    # but different from previous attempts
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 80.1285219342483
    
    rm(attempt_7)
    
## Attempt 8 -------------------------------------------------------------------
    # widen predation kernel even more, increase interaction to 1, 
    attempt_8 <- bonaire_model
    attempt_8 <- set_species_param_default(attempt_8, "sigma", 2.5)  
    attempt_8@interaction["predators", ] <- rep(1)
    attempt_8 <- reef_steady(attempt_8)
    
    # Converges immediately!
    attempt_8 <- calibrateReefBiomass(attempt_8)
    attempt_8 <- matchBiomasses(attempt_8)
    attempt_8 <- reef_steady(attempt_8)
    
    # Growth
    attempt_8 <- matchGrowth(attempt_8)
    attempt_8 <- reef_steady(attempt_8)
    
    # Growth - age mat still low for herbs
    age_mat_model = age_mat(attempt_8)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Iteration
    attempt_8 <- attempt_8 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Simulation run did not converge after 99 years. Value returned by the 
    # distance function was: 239.302546379538
    # Error in if (any(growth[idx] == 0)) { : 
    #         missing value where TRUE/FALSE needed
    
    # Distance same again - feels like something is causing it to oscillate
    rm(attempt_8)
    
## Attempt 9 -------------------------------------------------------------------
    # AHHHHHH I don't know what to try?
    # lower maturation age for predators?
    attempt_9 <- bonaire_model
    attempt_9 <- set_species_param_default(attempt_9, "sigma", 2.5)
    attempt_9@species_params$age_mat <- c(2,4,0.2)
    attempt_9@interaction["predators", ] <- rep(1)
    attempt_9 <- reef_steady(attempt_9)
    
    # Converges immediately!
    attempt_9 <- calibrateReefBiomass(attempt_9)
    attempt_9 <- matchBiomasses(attempt_9)
    attempt_9 <- reef_steady(attempt_9)
    
    # Growth
    attempt_9 <- matchGrowth(attempt_9)
    attempt_9 <- reef_steady(attempt_9)
    
    # Growth - age mat still low for herbs
    age_mat_model = age_mat(attempt_9)
    age_mat_observed = attempt_9@species_params$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Iteration
    attempt_9 <- attempt_9 |>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchGrowth()|> reef_steady()
    
    # Error in if (any(growth[idx] == 0)) { : 
    #         missing value where TRUE/FALSE needed
    
    # AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    
    