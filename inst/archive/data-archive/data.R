################################################################################
################################################################################

#' \linkS4class{MizerParams} object for a multi-species coral reef model
#'
#' Includes 10 species groups with biomass estimates from Karpata Reef
#' study site in the FORCE data set (Caribbean).
#'
#' PhD Thesis chapter 4
#' @keywords datasets
#'
#' @format A MizerParams object
#' @source PhD Thesis
"karpata_model"


#' species_params dataframe for a generic Caribbean reef
#'
#' Includes 10 species groups with biomass estimates from Karpata Reef
#' study site in the FORCE data set.
#'
#' @keywords datasets
#' PhD Thesis chapter 4
#'
#' @format data frame
#' @source PhD Thesis
"karpata_species"

#' Interaction matrix for a generic Caribbean reef
#'
#' Includes 10 species groups.
#'
#' @keywords datasets
#' PhD Thesis chapter 4
#'
#' @format data frame
#' @source PhD Thesis
"karpata_int"

#' Competitive method refuge parameters for a generic Caribbean reef
#'
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and `prop_protect` for each size bin, which gives the refuge density per
#' square meter for ten 5cm wide fish length bins ranging from 0 to 50 cm.
#' Data from the Karpata reef study site in Bonaire (FORCE data set).
#'
#' @keywords datasets
#' PhD Thesis chapter 4
#'
#' @format Data frame
#' @source PhD Thesis
"karpata_refuge"

#' Competitive method refuge parameters for a generic Caribbean reef
#'
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and `prop_protect` for each size bin, which gives the refuge density per
#' square meter for ten 5cm wide fish length bins ranging from 0 to 50 cm.
#' Data from the Aquarius reef study site.
#'
#' @keywords datasets
#' PhD Thesis chapter 4
#'
#' @format Data frame
#' @source PhD Thesis
"aquarius_refuge"

#' \linkS4class{MizerParams} object for a simple, generic Caribbean reef
#'
#' This object represents a simplified Caribbean reef ecosystem, parameterized
#' for use in mizerReef models. It includes three functional groups:
#' predators, herbivores, and invertebrates, with parameters derived from the
#' FORCE dataset (Future of Reefs in a Changing Environment,
#' see: http://theses.ncl.ac.uk/jspui/handle/10443/3574) and empirical
#' observations from Bonaire.
#'
#' The model structure and calibration are described in detail in Chapter 3 of
#' the PhD Thesis:
#' Beese, C. (2025). []"Modelling Coral Reef Futures: Exploring the role of
#' structural complexity in sustaining ecosystem services."](https://doi.org/10.26686/wgtn.26421523) Victoria University
#' of Wellington.
#'
#' @keywords datasets
#' @references FORCE dataset: http://theses.ncl.ac.uk/jspui/handle/10443/3574
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington. https://doi.org/10.26686/wgtn.26421523
#' @format A MizerParams object
#' @source PhD Thesis, FORCE dataset
"bonaire_model"

#' \linkS4class{MizerParams} object for a simple, generic Caribbean reef
#'
#' Includes 3 species groups: predators, herbivores, and invertebrates.
#'
#' @keywords datasets
#' PhD Thesis chapter 3 vignettes
#'
#' @format A MizerParams object
#' @source PhD Thesis
"bonaire_model"

#' species_params dataframe for a simple, generic Caribbean reef
#'
#' Includes 3 species groups: predators, herbivores, and invertebrates.
#'
#' @keywords datasets
#' PhD Thesis chapter 3 vignettes
#'
#' @format data frame
#' @source PhD Thesis
"bonaire_species"

#' interaction matrix for a simple, generic Caribbean reef
#'
#' Includes 3 species groups: predators, herbivores, and invertebrates.
#'
#' @keywords datasets
#' PhD Thesis chapter 3 vignettes
#'
#' @format data frame
#' @source PhD Thesis
"bonaire_interaction"

#' Competitive method refuge parameters for a simple, generic Caribbean reef
#'
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and `prop_protect` for each size bin, which gives the refuge density per
#' square meter for ten 5cm wide fish length bins ranging from 0 to 50 cm.
#'
#' @keywords datasets
#' PhD Thesis chapter 3 vignettes
#'
#' @format data frame
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
#' [newReefParams()] with the "binned" method and a proportional tuning profile.
#'
#' After species biomasses and growth rates have been calibrated to match
#' empirical observations, use the [newRefuge()] function to implement your
#' competitive refuge parameters. After using [newRefuge()], make sure to
#' iterate through [mizer:: matchBiomasses()], [matchReefGrowth()], and
#' [reefSteady()] again to regain the steady state.
#' @keywords datasets
#'
#' @format data frame
#' @source PhD Thesis
#' @concept calibration
"tuning_profile"

#' Stepped refuge profile for tuning steady states
#'
#' This is a 2-dimensional array containing start and end lengths for size bins
#' and `prop_protect` decreasing from 30% to 10% over the ten bins.
#'
#' This profile provides more protection to smaller size classes than larger
#' ones, as would be observed on a natural reef.
#'
#' These refuge parameters are intended for tuning the steady state when
#' using the density-dependent competitive method. The tuning profile provides
#' a constant proportion of refuges to all fish up to 50 cm in length.
#'
#' When creating a model using the competitive method, you should run
#' [newReefParams()] with the "binned" method and this tuning profile.
#'
#' After species biomasses and growth rates have been calibrated to match
#' empirical observations, use the [newRefuge()] function to implement your
#' competitive refuge parameters. After using [newRefuge()], make sure to
#' iterate through [mizer:: matchBiomasses()], [matchReefGrowth()], and
#' [reefSteady()] again to regain the steady state.
#' @keywords datasets
#'
#' @format data frame
#' @source PhD Thesis
#' @concept calibration
"step_tune"


#' Rubble trajectory refuge density scaling parameters
#'
#' This is a 2-dimensional array (refuge size x time post bleaching) containing
#' scaling values for refuge density for 15 years following bleaching
#'
#' @keywords datasets
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
#' @keywords datasets
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
#' @keywords datasets
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
#' @keywords datasets
"constant_scale"

#' Small example species parameters for body-shape/refuge method demo
#'
#' A small example species parameter data frame used in the Getting Started
#' vignette to demonstrate how different refuge methods affect vulnerability
#' profiles across species with different body shapes. The object contains one
#' row per functional group and includes columns used by `MizerReef` such as
#' `a`, `b`, `refuge_user`, `bad_pred`, and interaction preferences with
#' unstructured resources.
#'
#' @format A data frame with several rows and columns including:
#' \describe{
#'   \item{species}{character, group name}
#'   \item{l_max}{numeric, maximum length (cm)}
#'   \item{a,b}{numeric, length-weight parameters}
#'   \item{refuge_user}{logical, TRUE if group uses refuge}
#'   \item{bad_pred}{logical, TRUE if group accesses prey in refuge}
#'   \item{interaction_algae, interaction_detritus}{numeric 0/1 preference}
#' }
#' @source Vignette
#' @keywords datasets
#' @name body_shape_example_species_params
#' @docType data
"body_shape_example_species_params"

#' Small example interaction matrix for the body-shape demo
#'
#' Interaction matrix for the small example used in the Getting Started
#' vignette. Rows are predators and columns are prey groups. Values are
#' between 0 and 1 indicating interaction strength.
#'
#' @format A numeric matrix with row and column names
#' @source Vignette
#' @keywords datasets
#' @name body_shape_example_int
#' @docType data
"body_shape_example_int"
