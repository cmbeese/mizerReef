#' species_params dataframe for a simple, generic Caribbean reef
#'
#' This data frame contains species-level parameters for three functional
#' groups: predators, herbivores, and invertebrates. Parameters include growth
#' rates, mortality, and feeding traits, compiled as empirical data collected in
#' Bonaire from the FORCE dataset (see: http://theses.ncl.ac.uk/jspui/handle/10443/3574),
#' FishBase (www.fishbase.org), and published literature. See thesis for further
#' details on parameter sources and derivation.
#'
#' For details on parameter derivation and model calibration, see Chapter 3 of
#' the PhD Thesis:
#' Beese, C. (2025). ["Modelling Coral Reef Futures: Exploring the role of
#' structural complexity in sustaining ecosystem services."](https://doi.org/10.26686/wgtn.26421523) Victoria University
#' of Wellington.
#'
#' @keywords datasets
#' @references FORCE dataset: http://theses.ncl.ac.uk/jspui/handle/10443/3574
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington. https://doi.org/10.26686/wgtn.26421523
#' @format data frame
#' @source PhD Thesis, FORCE dataset
"bonaire_species"


#' interaction matrix for a simple, generic Caribbean reef
#'
#' This data frame encodes the trophic interactions among three functional
#' groups: predators, herbivores, and invertebrates. Interaction strengths are
#' parameterized using the FORCE dataset
#' (see: http://theses.ncl.ac.uk/jspui/handle/10443/3574) and empirical data
#' from Bonaire.
#'
#' For full details on interaction matrix construction, see Chapter 3 of the
#' PhD Thesis:
#' Beese, C. (2025). ["Modelling Coral Reef Futures: Exploring the role of
#' structural complexity in sustaining ecosystem services."](https://doi.org/10.26686/wgtn.26421523) Victoria University
#' of Wellington.
#'
#' @keywords datasets
#' @references FORCE dataset: http://theses.ncl.ac.uk/jspui/handle/10443/3574
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington. https://doi.org/10.26686/wgtn.26421523
#' @format data frame
#' @source PhD Thesis, FORCE dataset
"bonaire_interaction"


#' Competitive method refuge parameters for a simple, generic Caribbean reef
#'
#' This data frame contains refuge density parameters for a trait-based
#' competitive refuge model. It includes start and end lengths for size bins
#' and the proportion protected (`prop_protect`) for each bin, representing
#' refuge density per square meter for ten 5cm-wide fish length bins (0–50 cm).
#'
#' Parameters are derived from the FORCE dataset
#' (see: http://theses.ncl.ac.uk/jspui/handle/10443/3574) and field
#' measurements from Bonaire. For details on refuge parameterization
#' and model implementation, see Chapter 3 of the PhD Thesis:
#' Beese, C. (2025). ["Modelling Coral Reef Futures: Exploring the role of
#' structural complexity in sustaining ecosystem services."](https://doi.org/10.26686/wgtn.26421523) Victoria University
#' of Wellington.
#'
#' @keywords datasets
#' @references FORCE dataset: http://theses.ncl.ac.uk/jspui/handle/10443/3574
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington. https://doi.org/10.26686/wgtn.26421523
#' @format data frame
#' @source PhD Thesis, FORCE dataset
"bonaire_refuge"



#' Constant refuge profile for tuning steady states
#'
#' This data frame provides a constant refuge profile for model calibration and
#' steady-state tuning. It contains start and end lengths for size bins and sets
#' `prop_protect` to 20% for all bins up to 50 cm in length.
#'
#' These parameters are intended for use with the density-dependent competitive
#' method in mizerReef. For calibration, run [newReefParams()] with the "binned"
#' method and a proportional tuning profile. After calibrating species biomasses
#' and growth rates to match empirical observations (see FORCE dataset:
#' http://theses.ncl.ac.uk/jspui/handle/10443/3574), use [newRefuge()] to
#' implement competitive
#' refuge parameters. After updating, iterate through [mizer::matchBiomasses()],
#' [matchReefGrowth()], and [reefSteady()] to regain steady state.
#'
#' For full details on calibration and model setup, see Chapter 3 of the PhD
#' Thesis:
#' Beese, C. (2025). "Reef fish community dynamics and refuge availability: a
#' trait-based modeling approach." Victoria University of Wellington.
#' [https://doi.org/10.26686/thesis.123456](https://doi.org/10.26686/thesis.123456)
#'
#' @keywords datasets
#' @concept calibration
#' @references FORCE dataset: http://theses.ncl.ac.uk/jspui/handle/10443/3574
#' @references Beese, C. (2025). PhD Thesis. Victoria University of Wellington. [https://doi.org/10.26686/thesis.123456](https://doi.org/10.26686/thesis.123456)
#' @format data frame
#' @source PhD Thesis, FORCE dataset
"tuning_profile"
