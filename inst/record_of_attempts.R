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

#bonaire_species$interaction_resource <- c(1,0,1)
bonaire_species$gamma <- c(6.4,0.2,0.2)

# Set model --------------------------------------------------------------------
bonaire_model <- newReefParams(group_params = bonaire_species,
                               interaction = bonaire_int,
                               method = "competitive",
                               method_params = bonaire_refuge)
## Attempt 4 -------------------------------------------------------------------
    # widening the predation kernel, setting all interactions to 1 and 
    # setting gamma to value from Alice's model - 6.4 for preds, 0.2 inv and herb
    attempt_4 <- bonaire_model
    # attempt_4 <- set_species_param_default(attempt_4, "sigma", 2)
    #attempt_4@species_params$gamma <- c(6.4,0.2,0.2)
    attempt_4@interaction["predators", ] <- rep(1)
    attempt_4 <- reef_steady(attempt_4)
    
    # Converges immediately!
    attempt_4 <- calibrateReefBiomass(attempt_4)
    attempt_4 <- matchBiomasses(attempt_4)
    attempt_4 <- reef_steady(attempt_4)
    
    # Growth
    attempt_4 <- matchReefGrowth(attempt_4)
    attempt_4 <- reef_steady(attempt_4)
    
    # Iteration
    attempt_4 <- attempt_4 |>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()
    
    # Maturation age 
    age_mat_model = age_mat(attempt_4)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    # Look at spectra - look crazy!
    plotSpectra(attempt_4, power = 2)
    
    # Resource is low - increase resource
    attempt_4 <- scaleDownBackground(attempt_4, factor = 1/2)
    
    # Find steady again
    attempt_4 <- attempt_4 |> matchReefGrowth() |> reef_steady()
    
    # Maturation age 
    age_mat_model = age_mat(attempt_4)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    plotSpectra(attempt_4, power = 2)
    
    # Increase resource again
    attempt_4 <- scaleDownBackground(attempt_4, factor = 1/2)
    
    # Find steady again
    attempt_4 <- attempt_4 |> matchReefGrowth() |> reef_steady()
    plotSpectra(attempt_4, power = 2)
    
    # Reproduction
    attempt_4 <- setBevertonHolt(attempt_4, erepro = 0.03)
    
    attempt_4 <- attempt_4 |>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()
    
    plotSpectra(attempt_4, power = 2)
    # Herbivores seem to be going crazy?
    
    # Increase resource again
    attempt_4 <- scaleDownBackground(attempt_4, factor = 1/2)
    # Find steady again
    attempt_4 <- attempt_4 |>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()
    
    # Plot spectra again
    plotSpectra(attempt_4, power = 2)

    rm(attempt_4)
    
## Attempt 6 -------------------------------------------------------------------
# widening the predation kernel, setting all interactions to 1, and
# increasing gamma to higher values from range in Nash 2015, let invertebrates 
# get bigger so predators can feed on them as well
    attempt_6 <- bonaire_model
    #attempt_6 <- set_species_param_default(attempt_6, "sigma", 2)  
    attempt_6@interaction["predators", ] <- rep(1)
    attempt_6@species_params$l_max <- c(50,50,30)
    attempt_6@species_params$gamma <- c(6.4,0.2,0.2)
    #attempt_6@species_params$gamma <- c(100,NA,NA)
    attempt_6 <- reef_steady(attempt_6)
    
    
    p<- attempt_6@species_params
    
    # Converges immediately!
    attempt_6 <- calibrateReefBiomass(attempt_6)
    attempt_6 <- matchBiomasses(attempt_6)
    attempt_6 <- reef_steady(attempt_6)
    
    # Growth
    attempt_6 <- matchReefGrowth(attempt_6)
    attempt_6 <- reef_steady(attempt_6)
    
    # Iteration
    attempt_6 <- attempt_6 |>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> 
        matchReefGrowth()|> reef_steady()
    
    # Resource looks low - try raising it match sheldon's spectrum
    plotSpectra(attempt_6, total = TRUE, power = 1)
    plotSpectra(attempt_6, total = TRUE, power = 2)
    plotVulnerable(attempt_6)
    
    attempt_7 <- scaleDownBackground(attempt_6, factor = 1/2)
    
    attempt_7 <- attempt_7 |>
        matchBiomasses()|> matchReefGrowth()|> reef_steady()|>
        matchBiomasses()|> matchReefGrowth()|> reef_steady()|>
        matchBiomasses()|> matchReefGrowth()|> reef_steady()|>
        matchBiomasses()|> matchReefGrowth()|> reef_steady()|>
        matchBiomasses()|> matchReefGrowth()|> reef_steady()|>
        matchBiomasses()|> matchReefGrowth()|> reef_steady()

    plotSpectra(attempt_7, total = TRUE, power = 2)
    
    
  
    
    rm(attempt_7)
    plotDiet(attempt_7)
    plotFeedingLevel(attempt_6)
    plotDiet(attempt_6)
    # Maturation age for herbivores still very low
    age_mat_model = age_mat(attempt_6)
    age_mat_observed = bonaire_species$age_mat
    data.frame(age_mat_model, age_mat_observed)
    
    rep_level <- getReproductionLevel(attempt_6)

    
    attempt_6
    
    rm(attempt_6)
    
    