#' species_params dataframe for a simple, generic Caribbean reef with 3
#' species groups
#'
#' This data frame contains species-level parameters for three functional
#' groups: predators, herbivores, and invertebrates. Parameters include
#' biomasses from empirical data collected in Bonaire (FORCE dataset,
#' see: Dryden, C. (2017), "Habitat structural complexity of Caribbean
#' coral reefs and its relationships with fish community structure",
#' Newcastle University), as well as growth
#' rates, maturation ages, and other life-history traits derived from FishBase
#' (www.fishbase.org) and published literature. See thesis for further details
#' on parameter sources and derivation.
#'
#' For details on parameter derivation and model calibration, see Chapter 3 of
#' the PhD Thesis:
#' Beese, C. (2025). "Modelling Coral Reef Futures: Exploring the role of
#' structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington.
#'
#' @keywords datasets
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' http://theses.ncl.ac.uk/jspui/handle/10443/3574
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem
#' services." Victoria University of Wellington.
#' https://doi.org/10.26686/wgtn.26421523
#' @format data frame
#' @source PhD Thesis, FORCE dataset
"caribbean_3_species"


#' interaction matrix for a simple, generic Caribbean reef with 3 species
#' groups
#'
#' This data frame encodes the trophic interactions among three functional
#' groups: predators, herbivores, and invertebrates. The interaction
#' strengths set to 1 for all piscivores in this model. Herbivores consume
#' algae and invertebrates consume both plankton and detritus, with the
#' proportion of detritus increasing with size.
#'
#' For full details on interaction matrix construction, see Chapter 3 of the
#' PhD Thesis:
#' Beese, C. (2025). "Modelling Coral Reef Futures: Exploring the role of
#' structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington.
#'
#' @keywords datasets
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem
#' services." Victoria University of Wellington.
#' https://doi.org/10.26686/wgtn.26421523
#' @format data frame
#' @source PhD Thesis, FORCE dataset
"caribbean_3_interaction"


#' mizerParams object for a simple, generic Caribbean reef with 3 species groups
#'
#' A `MizerParams` object for a simple
#' Caribbean reef ecosystem model with three functional groups: predators,
#' herbivores, and invertebrates. It is built using the `caribbean_3_species`
#' and `caribbean_3_interaction` data, and includes all parameter settings and
#' calibrations as described in the associated PhD thesis and FORCE dataset.
#' This object can be used directly with mizer functions for simulation and
#' analysis.
#'
#' For details on parameter derivation, model calibration, and implementation,
#' see Chapter 3 of the PhD Thesis:
#' Beese, C. (2025). "Modelling Coral Reef Futures: Exploring the role of
#' structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington.
#'
#' @keywords datasets
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish community
#' structure. PhD Thesis, Newcastle University.
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington.
#' https://doi.org/10.26686/wgtn.26421523
#' @format mizerParams object
#' @source PhD Thesis, FORCE dataset
"caribbean_3_model"


#' Competitive method refuge parameters for Karpata reef in Bonaire
#'
#' This data frame contains refuge density parameters for a trait-based
#' competitive refuge model. It includes start and end lengths for size bins
#' and the proportion protected (`prop_protect`) for each bin, representing
#' refuge density per square meter for ten 5cm-wide fish length bins
#' (0–50 cm).
#'
#' Parameters are derived from the FORCE dataset
#' (see: Dryden, C. (2017)) and field
#' measurements from Bonaire. For details on refuge parameterization and
#' model implementation, see Chapter 3 of the PhD Thesis:
#' Beese, C. (2025). "Modelling Coral Reef Futures: Exploring the role
#' of structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington.
#'
#' @keywords datasets
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem
#' services." Victoria University of Wellington.
#' https://doi.org/10.26686/wgtn.26421523
#' @format data frame
#' @source PhD Thesis, FORCE dataset
"karpata_refuge"


#' Constant refuge profile for tuning steady states
#'
#' This data frame provides a constant refuge profile for model calibration
#' and steady-state tuning. It contains start and end lengths for size bins
#' and sets `prop_protect` to 60% for all bins up to 50 cm in length.
#'
#' These parameters are intended to be used with the "binned" method, which
#' is independent of density. The tuning profile should be used while
#' calibrating biomass and growth rates before switching to the competitive
#' method. See the "Getting Started" vignette for more details.
#'
#' @keywords datasets
#' @concept calibration
# @seealso \link{Getting Started}
#' @references Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
#' [https://doi.org/10.26686/thesis.123456](https://doi.org/10.26686/thesis.123456)
#' @format data frame
#' @source Beese PhD Thesis
"tuning_profile"
