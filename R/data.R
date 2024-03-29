# \linkS4class{MizerParams} object for testing
# 
# PhD Thesis chapter 4
#
# @format dataframe
# @source PhD Thesis
#"test"

# species_params dataframe for testing
# 
# PhD Thesis chapter 4
#
# @format dataframe
# @source PhD Thesis
# "test_sp"

# interaction dataframe for testing
# 
# PhD Thesis chapter 4
#
# @format dataframe
# @source PhD Thesis
# "test_i"

################################################################################
################################################################################

#' \linkS4class{MizerParams} object for a multispecies generic reef
#' 
#' Includes 10 species groups with biomass estimates from Karpata Reef 
#' study site in the FORCE data set (Caribbean).
#' 
#' PhD Thesis chapter 4
#'
#' @format A MizerParams object
#' @source PhD Thesis
"karpata_model"


#' species_params dataframe for a generic Caribbean reef
#' 
#' Includes 10 species groups with biomass estimates from Karpata Reef 
#' study site in the FORCE data set.
#' 
#' PhD Thesis chapter 4
#'
#' @format dataframe
#' @source PhD Thesis
"karpata_species"

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
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and `prop_protect` for each size bin, which gives the refuge density per 
#' square meter for ten 5cm wide fish length bins ranging from 0 to 50 cm.
#' Data from the Karpata reef study site.
#'  
#' PhD Thesis chapter 4
#'
#' @format A MizerParams object
#' @source PhD Thesis
"karpata_refuge"

#' Competitive method refuge parameters for a generic Caribbean reef
#' 
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and `prop_protect` for each size bin, which gives the refuge density per 
#' square meter for ten 5cm wide fish length bins ranging from 0 to 50 cm.
#' Data from the Aquarius reef study site.
#' 
#' PhD Thesis chapter 4
#'
#' @format A MizerParams object
#' @source PhD Thesis
"aquarius_refuge"

#' \linkS4class{MizerParams} object for a simple, generic Caribbean reef
#' 
#' Includes 3 species groups: predators, herbivores, and invertebrates. 
#' 
#' PhD Thesis chapter 3 vignettes
#'
#' @format A MizerParams object
#' @source PhD Thesis
"bonaire_model"

#' species_params dataframe for a simple, generic Caribbean reef
#' 
#' Includes 3 species groups: predators, herbivores, and invertebrates. 
#' 
#' PhD Thesis chapter 3 vignettes
#' 
#' @format dataframe
#' @source PhD Thesis
"bonaire_species"

#' interaction matrix for for a simple, generic Caribbean reef
#' 
#' Includes 3 species groups: predators, herbivores, and invertebrates. 
#' 
#' PhD Thesis chapter 3 vignettes
#' 
#' @format dataframe
#' @source PhD Thesis
"bonaire_int"

#' Competitive method refuge parameters for a simple, generic Caribbean reef
#' 
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and `prop_protect` for each size bin, which gives the refuge density per 
#' square meter for ten 5cm wide fish length bins ranging from 0 to 50 cm.
#' 
#' PhD Thesis chapter 3 vignettes
#' 
#' @format dataframe
#' @source PhD Thesis
"bonaire_refuge"

#' Constant refuge profile for tuning steady states
#' 
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and `prop_protect` equal to 20% for all size bins up to 50 cm in length.
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
#' @concept calibration
"tuning_profile"

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
#' @concept calibration
"step_tune"


#' Rubble trajectory refuge density scaling parameters
#' 
#' This is a 2-dimensional array (refuge size x time post bleaching) containing 
#' scaling values for refuge density for 15 years following bleaching
#' 
#' PhD Thesis chapter 5
#' 
#' @concept degradation
#' @format matrix
#' @source PhD Thesis
"rubble_scale"

#' Algae trajectory refuge density scaling parameters
#' 
#' This is a 2-dimensional array (refuge size x time post bleaching) containing 
#' scaling values for refuge density for 15 years following bleaching
#' 
#' PhD Thesis chapter 5
#' 
#' @concept degradation
#' @format matrix
#' @source PhD Thesis
"algae_scale"

#' Recovery trajectory refuge density scaling parameters
#' 
#' This is a 2-dimensional array (refuge size x time post bleaching) containing 
#' scaling values for refuge density for 15 years following bleaching
#' 
#' PhD Thesis chapter 5
#' 
#' @concept degradation
#' @format matrix
#' @source PhD Thesis
"recovery_scale"

#' Trajectory with no refuge density scaling for testing
#' 
#' This is a 2-dimensional array (refuge size x time post bleaching) containing 
#' scaling values for refuge density for 15 years with no degradation. 
#' For testing.
#' 
#' PhD Thesis chapter 5
#' 
#' @concept degradation
#' @format matrix
#' @source PhD Thesis
"constant_scale"

