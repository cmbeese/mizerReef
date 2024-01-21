# Setting up a Caribbean coral reef mizer model with multiple resources
# Model steady state calibration
# last tuned 21st January 2024

## Setup - load packages -------------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)
library(here)

## Load parameters -------------------------------------------------------------

    # Load species parameter data
    bonaire_species <- read.csv(here("inst/bonaire_species.csv"))
    bonaire_int     <- read.csv(here("inst/bonaire_int.csv"),  row.names = 1)
    bonaire_refuge  <- bonaire_refuge
    
    # With these parameters, herbivores consume plankton at small sizes and 
    #   transition fully to algae by maturity 
    # With these parameters, invertebrates consume plankton and detritus,
    #   with the proportion of detritus increasing with size

## Set model -------------------------------------------------------------------
    params <- newReefParams(group_params = bonaire_species,
                            interaction = bonaire_int,
                            method = "competitive",
                            method_params = bonaire_refuge)
    
## Project to first steady state -----------------------------------------------
    params <- reef_steady(params)

## Calibrate biomasses and growth ----------------------------------------------
    
    # Match observed species group biomasses
    params <- calibrateReefBiomass(params)
    params <- matchBiomasses(params)
    params <- reef_steady(params)
    
    # Match observed growth rates
    params <- matchReefGrowth(params)
    params <- reef_steady(params)
    
    # Check for match with age at maturity
    age_mat_observed = bonaire_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)
    # These are already quite close 
    
    # Check biomass match - still way off
    plotBiomassVsSpecies(params)
    
## Tune home range -------------------------------------------------------------
    
    # Check current search volume
    current_vol <- getSearchVolume(params)
    g <- params@species_params$gamma
    g[1] <- 1
    # By Nash 2014, gamma for predators should be at least 1 for predators and
    # at least 0.004 for herbivores
    # Looks okay for herbs but should be higher for preds
    
    # Implement new search rate
    params@species_params$gamma <- g
    gp <- 1
    q <- 0.75
    pred_search <- gp*(params@w^q)
    herb_search <- current_vol["herbivores",]
    inv_search  <- current_vol["inverts",]
    new_vol <- rbind(pred_search, herb_search, inv_search)
    rownames(new_vol) <- c("predators","herbivores","inverts")
    params <- setSearchVolume(params, search_vol = new_vol)
    
## Steady state iteration ------------------------------------------------------
    
    # Iterate to refine growth and biomass
    params <- params |>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()
    
    # Much better
    plotBiomassVsSpecies(params)
    
    # Check match with observed age at maturity
    age_mat_observed = bonaire_species$age_mat
    age_mat_model = age_mat(params)
    data.frame(age_mat_model, age_mat_observed)
    # Closer than needed
    
## Check resulting spectra and tune resources ----------------------------------
    
    # Resource looks low - should match sheldon's spectrum
    # looks fairly straight not bad but some bumps
    plotSpectra(params, total = TRUE, power = 1)
    plotSpectra(params, total = TRUE, power = 2)
    
    # plot feeding level to check if resource is too low
    plotFeedingLevel(params, species = "inverts")
    
    # Invert feeding level is relatively stable through life, non-linearities
    #   are probably due to refuge
    
# Tune reproduction ------------------------------------------------------------
    # We do not have yield or catch data - can't tune size distribution
    # First attempt to set very low to see what the minimum values are
    params <- setBevertonHolt(params, erepro = 0.0001)
    # Now set setting erepro same for all species, as low as possible
    params <- setBevertonHolt(params, erepro = 0.38)
    # Check reproduction level (value between 0 and 1) - should be higher for
    # larger, slow growing species and low for small, fast growing ones
    rep <- getReproductionLevel(params)
    # These are low for predators and herbivores, and strangely high for 
    # inverts. A reproduction level closer to one means reproduction rate is 
    # almost totally independent of the investment into reproduction
    # Reproduction should be density independent on reefs
    
    # Check comparison of density dependent & independent reduction
    getRDI(params) / getRDD(params)
    # Reproduction is density independent for inverts, density dependent for 
    # preds and herbs - but possible too much - RDI is only slightly higher
    
    # Let's increase reproduction level to 0.5 for predators and herbivores
    # so that 
    rep_level <- c(0.5, 0.5, rep[3])
    names(rep_level) <- c("predators","herbivores","inverts")
    params <- setBevertonHolt(params,
                              reproduction_level = rep_level)
    
    # Iterate to get back to steady state
    params <- params |>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()|>
        calibrateReefBiomass()|> matchBiomasses()|> matchReefGrowth()|> 
        reef_steady()
    
    # Check new reproduction - these look better
    rep <- getReproductionLevel(params)
    getRDI(params) / getRDD(params)
    
    # Check new spectra
    plotSpectra(params, total = TRUE, power = 1)
    plotSpectra(params, total = TRUE, power = 2)
    
    # Save!
    bonaire_model <- reef_steady(params)
    
# Plots ------------------------------------------------------------------------
plotBiomassVsSpecies(bonaire_model)
plotRefuge(bonaire_model)
plotSpectra(bonaire_model, power = 1, total = TRUE)
plotDiet(bonaire_model)  

# I am happy with these parameters!
# Save object ------------------------------------------------------------------

    # Params object
    save(bonaire_model,   file = "data/bonaire_model.rda")

    # CSV Files
    save(bonaire_species, file = "data/bonaire_species.rda")
    save(bonaire_int,     file = "data/bonaire_int.rda")
    save(bonaire_refuge,  file = "data/bonaire_refuge.rda")
    