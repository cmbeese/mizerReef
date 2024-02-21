#' species_params dataframe for a generic Caribbean reef
#' 
#' Biomass estimates from Karpata Reef 
#' 
#' Includes 10 species groups.
#' 
#' PhD Thesis chapter 4
#'
#' @format dataframe
#' @source PhD Thesis
"karpata_species"

# species_params dataframe for a generic Caribbean reef
# 
# Biomass estimates from Aquarius reef.
# 
# Includes 10 species groups.
# 
# PhD Thesis chapter 4
#
# @format dataframe
# @source PhD Thesis
#"aquarius_species"

#' Interaction matrix for a generic Caribbean reef
#' 
#' Includes 10 species groups.
#' 
#' PhD Thesis chapter 4
#'
#' @format dataframe
#' @source PhD Thesis
"karpata_int"

#' Competitive method refuge parameters for a generic Caribbean reef
#' 
#' Refuge densities from Karpata Reef site
#'  
#' PhD Thesis chapter 4
#'
#' @format A MizerParams object
#' @source PhD Thesis
"karpata_refuge"

#' Competitive method refuge parameters for a generic Caribbean reef
#' 
#' Refuge densities from Aquarius Reef site
#' 
#' PhD Thesis chapter 4
#'
#' @format A MizerParams object
#' @source PhD Thesis
"aquarius_refuge"

#' MizerParams object for a simple, generic Caribbean reef
#' 
#' Includes predators, herbivores, and invertebrates. 
#' 
#' PhD Thesis chapter 3 vignettes
#'
#' @format A MizerParams object
#' @source PhD Thesis
"bonaire_model"

#' species_params for a simple, generic Caribbean reef
#' 
#' Includes predators, herbivores, and invertebrates. 
#' 
#' PhD Thesis chapter 3 vignettes
#' 
#' @format dataframe
#' @source PhD Thesis
"bonaire_species"

#' interaction matrix for for a simple, generic Caribbean reef
#' 
#' Includes predators, herbivores, and invertebrates. 
#' 
#' PhD Thesis chapter 3 vignettes
#' 
#' @format dataframe
#' @source PhD Thesis
"bonaire_int"

#' refuge parameters for for a simple, generic Caribbean reef
#' 
#' Includes predators, herbivores, and invertebrates. 
#' 
#' PhD Thesis chapter 3 vignettes
#' 
#' @format datafram
#' @source PhD Thesis
"bonaire_refuge"

#' Constant refuge profile for tuning steady states
#' 
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and prop_protect equal to 20% for all size bins up to 50 cm in length.
#' 
#' These refuge parameters are intended for tuning the steady state when
#' using the density-dependent competitive method. The tuning profile provides
#' a constant proportion of refuges to all fish up to 50 cm in length. 
#' 
#' When creating a model using the competitive method, you should run
#' [newReefParams()] with the "binned" method and a proportional 
#' tuning profile.
#' 
#' After species biomasses and growth rates have been calibrated to match
#' empirical observations, use the [newRefuge()] function to implement your
#' competitive refuge parameters. After using [newRefuge()], make sure to
#' iterate through [mizer:: matchBiomasses()], [matchReefGrowth()], and 
#' [reefSteady()] again to regain the steady state.
#' 
#' @format dataframe         
#' @source PhD Thesis
"constant_tune"

#' Stepped refuge profile for tuning steady states
#' 
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and prop_protect decreasing from 30% to 10% over the ten bins. 
#' 
#' This profile provides more protection to smaller size classes than larger
#' ones, as would be observed on a natural reef. 
#' 
#' These refuge parameters are intended for tuning the steady state when
#' using the density-dependent competitive method. The tuning profile provides
#' a constant proportion of refuges to all fish up to 50 cm in length. 
#' 
#' When creating a model using the competitive method, you should run
#' [newReefParams()] with the "binned" method and thing tuning profile.
#' 
#' After species biomasses and growth rates have been calibrated to match
#' empirical observations, use the [newRefuge()] function to implement your
#' competitive refuge parameters. After using [newRefuge()],  make sure to
#' iterate through [mizer:: matchBiomasses()], [matchReefGrowth()], and 
#' [reefSteady()] again to regain the steady state.
#' 
#' @format dataframe         
#' @source PhD Thesis
"step_tune"