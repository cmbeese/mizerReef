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
#' Two herbivore parameters were revised in 2026 during recalibration against
#' a corrected senescence-mortality formula (see `caribbean_3_model`'s
#' documentation and `inst/scripts/Caribbean_3_model-calibration.R`'s design
#' note): `satiation` (FALSE -> TRUE) and `age_mat` (4 -> 1.6, the age at
#' median sexual maturity reported for the stoplight parrotfish, *Sparisoma
#' viride*, by Rivera Hernandez & Shervette (2025), replacing a value that
#' appears to have conflated it with that species' age at sexual transition).
#'
#' @keywords datasets
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' http://theses.ncl.ac.uk/jspui/handle/10443/3574
#' @references Rivera Hernandez, J.M. & Shervette, V.R. (2025). Addressing
#' life history information gaps for Caribbean parrotfishes: queen parrotfish
#' Scarus vetula and stoplight parrotfish Sparisoma viride. Environmental
#' Biology of Fishes, 108, 179-198. https://doi.org/10.1007/s10641-024-01651-x
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
#' Recalibrated in 2026 against a corrected senescence-mortality formula
#' (`reefSenMort()`) -- fish abundances, growth/consumption rates and
#' reproduction parameters were all re-tuned to restore a genuine steady
#' state under the corrected code (verified: repeated `reefSteady()` calls
#' leave biomass unchanged). Two herbivore species_params were also changed
#' from their original thesis values as part of this: `satiation` is now
#' `TRUE` (was `FALSE`) since herbivore biomass has no density-dependent
#' brake without some cap on individual intake once mortality is realistically
#' low -- the realised feeding level is consistently close to 1, still
#' consistent with the citations behind the package's usual herbivore
#' default (Caribbean herbivores' guts observed to be full nearly
#' continuously); and `age_mat` is now `1.6` (was `4`), the age at median
#' sexual maturity for *Sparisoma viride* (Rivera Hernandez & Shervette
#' 2025), replacing a value that appears to have conflated it with that
#' species' age at sexual transition. See
#' `inst/scripts/Caribbean_3_model-calibration.R`'s design note and
#' `vignettes/caribbean_3_model-description.Rmd` for the full reasoning.
#' Biomasses (predators/herbivores/inverts, g/m^2) are now 107.9/34.1/40.1,
#' matching the FORCE-survey targets (107/34/40) essentially exactly.
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
#' This data frame contains refuge density parameters for the competitive
#' refuge method. It records the density of refuges (no./m^2) in each of
#' ten 5 cm fish-length bins (0-50 cm), derived from field measurements at
#' the Karpata reef site in Bonaire (FORCE dataset, Dryden 2017).
#'
#' Use this with \code{method = "competitive"} in \code{\link{setRefuge}()}.
#'
#' This profile is reused across mizerReef's example models: it is the
#' data behind \code{\link{caribbean_10_model}}'s refuge as well as
#' \code{\link{caribbean_3_model}}'s.
#'
#' Use this with \code{method = "competitive"} in \code{\link{setRefuge}()}.
#'
#' @keywords datasets
#' @concept refugeParams
#' @seealso \code{\link{setRefuge}()}, \code{\link{caribbean_10_model}},
#'   \code{\link{caribbean_3_model}}, \code{\link{aquarius_refuge}}
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem
#' services." Victoria University of Wellington.
#' https://doi.org/10.26686/wgtn.26421523
#' @format A data frame with 10 rows and 3 columns: \code{start_L} (start of
#'   length bin in cm), \code{end_L} (end of length bin in cm), and
#'   \code{refuge_density} (refuge density in no./m^2).
#' @source FORCE dataset, Bonaire (Dryden 2017)
"karpata_refuge"


#' Competitive method refuge parameters for Aquarius reef research station
#'
#' This data frame contains refuge density parameters for the competitive
#' refuge method. It records the density of refuges (no./m^2) in each of
#' ten 5 cm fish-length bins (0-50 cm), derived from field measurements at
#' the Aquarius reef research station (FORCE dataset, Dryden 2017), from
#' the same field hole-density survey pipeline as \code{\link{karpata_refuge}}.
#' Like \code{karpata_refuge}, this profile is general-purpose and can be
#' used with any mizerReef model, not just one particular example model.
#'
#' Use this with \code{method = "competitive"} in \code{\link{setRefuge}()}.
#'
#' @keywords datasets
#' @concept refugeParams
#' @seealso \code{\link{setRefuge}()}, \code{\link{karpata_refuge}}
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' @references Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures:
#' Exploring the role of structural complexity in sustaining ecosystem
#' services." Victoria University of Wellington.
#' https://doi.org/10.26686/wgtn.26421523
#' @format A data frame with 10 rows and 3 columns: \code{start_L} (start of
#'   length bin in cm), \code{end_L} (end of length bin in cm), and
#'   \code{refuge_density} (refuge density in no./m^2).
#' @source FORCE dataset, Aquarius reef research station (Dryden 2017)
"aquarius_refuge"


#' Species parameters for a generic Caribbean reef (10 functional groups)
#'
#' A species parameter data frame for a 10-group size-spectrum model of
#' a Caribbean coral reef. Functional groups include multiple predator
#' guilds (engulf, grab, eel-like, cryptic, invertebrate, planktivorous),
#' parrotfish, farming damselfish, other herbivores, and invertebrates.
#' Parameters are derived from field measurements at Karpata Reef, Bonaire
#' (FORCE dataset, Dryden 2017) and published life-history data (FishBase).
#'
#' @keywords datasets
#' @seealso \code{\link{caribbean_10_model}}, \code{\link{caribbean_10_interaction}},
#'   \code{\link{karpata_refuge}}
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' @references Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
#'   https://doi.org/10.26686/wgtn.26421523
#' @format A data frame with one row per functional group and columns for
#'   life-history parameters (w_max, w_mat, age_mat, h, beta, sigma, etc.),
#'   refuge use (\code{refuge_user}, \code{blocked_pred}, \code{satiation}),
#'   and resource interactions (\code{interaction_algae},
#'   \code{interaction_detritus}).
#' @source PhD Thesis, FORCE dataset
"caribbean_10_species"


#' Interaction matrix for a generic Caribbean reef (10 functional groups)
#'
#' The predator-prey interaction matrix for the 10-group Caribbean reef model.
#' Rows represent predators, columns prey. Values encode the relative
#' vulnerability of each prey group to each predator.
#'
#' @keywords datasets
#' @seealso \code{\link{caribbean_10_model}}, \code{\link{caribbean_10_species}}
#' @references Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
#'   https://doi.org/10.26686/wgtn.26421523
#' @format A numeric matrix (10 x 10).
#' @source PhD Thesis, FORCE dataset
"caribbean_10_interaction"


#' mizerReef model for a generic Caribbean reef (10 functional groups)
#'
#' A calibrated \code{mizerReef} object for a generic Caribbean coral reef
#' with ten functional groups, parameterised from field data collected at
#' Karpata Reef, Bonaire (FORCE dataset, Dryden 2017). Functional groups span
#' six predator guilds, parrotfish, farming damselfish, other herbivores, and
#' invertebrates. The model includes competitive predation refuge (using
#' \code{\link{karpata_refuge}}) and benthic resources (algae and detritus).
#'
#' Built from \code{\link{caribbean_10_species}},
#' \code{\link{caribbean_10_interaction}}, and \code{\link{karpata_refuge}};
#' calibrated to match observed biomasses and upgraded to mizerReef 2.1.0 /
#' mizer 3.1.0. Detritus resource dynamics are tuned via [tuneUR()] so that
#' they are genuinely at steady state (\code{dB/dt = 0}) for the model's
#' calibrated abundances. Algae biomass is likewise tuned by [tuneUR()] to
#' its own steady state, but for a fixed, literature-informed production
#' rate (\code{algae_growth}, left at \code{newReefParams()}'s default of
#' \code{2000} grams per square meter per year -- see
#' \code{\link{setAlgaeParams}}) that is not itself retuned to match
#' consumption, since real algal production is not driven by grazer demand;
#' an earlier untuned version of this object (algae biomass not yet solved
#' for its steady state) is archived at
#' \code{inst/archive/caribbean_10_model_untuned.rda} for reference.
#'
#' @keywords datasets
#' @seealso \code{\link{caribbean_10_species}}, \code{\link{caribbean_10_interaction}},
#'   \code{\link{karpata_refuge}}, \code{\link{newReefParams}}
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' @references Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
#'   https://doi.org/10.26686/wgtn.26421523
#' @format A \code{mizerReef} object (extends [MizerParams]) with
#'   10 species, competitive refuge, and algae/detritus components.
#' @source PhD Thesis, FORCE dataset
"caribbean_10_model"


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


#' Degradation scaling matrix for the rubble trajectory
#'
#' A 10 x 15 numeric matrix of multiplicative scaling factors for refuge
#' density under the rubble disturbance trajectory. Each row corresponds to
#' a 5 cm refuge size bin (0-5 cm through 45-50 cm) and each column to a
#' simulation year (column 1 = bleaching year, columns 2-15 = post-bleaching
#' years 1-14). Values are applied as \code{new_rd = scale * old_rd} inside
#' \code{\link{reefDegrade}()}, so a value of 1.4 means a 40\% increase in
#' refuge density and 0.4 means a 60\% reduction.
#'
#' In the rubble trajectory, large structural refuges collapse while an
#' initial pulse of rubble temporarily increases the availability of small
#' refuges (0-5 cm bin) before all sizes return to baseline.
#'
#' Values are derived from \code{inst/data-csv/deg_scales.csv} by computing
#' \code{multiplier = 1 + delta} for each cell.
#'
#' @keywords datasets
#' @concept degradation
#' @seealso \code{\link{algae_scale}}, \code{\link{recovery_scale}},
#'   \code{\link{setDegradation}}, \code{\link{reefDegrade}}
#' @format A numeric matrix with 10 rows (refuge size bins) and 15 columns
#'   (time steps). Row names are size bin labels (e.g. \code{"0 to 5"});
#'   column names are integers 1 to 15.
#' @references Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
#'   https://doi.org/10.26686/wgtn.26421523
"rubble_scale"


#' Degradation scaling matrix for the algae trajectory
#'
#' A 10 x 15 numeric matrix of multiplicative scaling factors for refuge
#' density under the algae disturbance trajectory. Each row corresponds to
#' a 5 cm refuge size bin (0-5 cm through 45-50 cm) and each column to a
#' simulation year (column 1 = bleaching year, columns 2-15 = post-bleaching
#' years 1-14).
#'
#' In the algae trajectory, moderate coral loss is partially offset by algal
#' colonisation of rubble. Small refuges (0-5 cm bin) show a slight
#' overshoot above baseline (1.02x) during years 6-9 before returning to
#' pre-disturbance levels.
#'
#' Values are derived from \code{inst/data-csv/deg_scales.csv} by computing
#' \code{multiplier = 1 + delta} for each cell.
#'
#' @keywords datasets
#' @concept degradation
#' @seealso \code{\link{rubble_scale}}, \code{\link{recovery_scale}},
#'   \code{\link{setDegradation}}, \code{\link{reefDegrade}}
#' @format A numeric matrix with 10 rows (refuge size bins) and 15 columns
#'   (time steps). Row names are size bin labels; column names are integers
#'   1 to 15.
#' @references Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
#'   https://doi.org/10.26686/wgtn.26421523
"algae_scale"


#' Degradation scaling matrix for the recovery trajectory
#'
#' A 10 x 15 numeric matrix of multiplicative scaling factors for refuge
#' density under a recovery trajectory. Each row corresponds to a 5 cm
#' refuge size bin (0-5 cm through 45-50 cm) and each column to a simulation
#' year (column 1 = bleaching year, columns 2-15 = post-bleaching years
#' 1-14).
#'
#' In the recovery trajectory, reefs experience an initial decline in refuge
#' availability followed by a sustained partial recovery above pre-disturbance
#' levels. Small and medium bins recover more quickly than large bins.
#'
#' Values are derived from \code{inst/data-csv/deg_scales.csv} by computing
#' \code{multiplier = 1 + delta} for each cell.
#'
#' @keywords datasets
#' @concept degradation
#' @seealso \code{\link{rubble_scale}}, \code{\link{algae_scale}},
#'   \code{\link{setDegradation}}, \code{\link{reefDegrade}}
#' @format A numeric matrix with 10 rows (refuge size bins) and 15 columns
#'   (time steps). Row names are size bin labels; column names are integers
#'   1 to 15.
#' @references Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
#'   https://doi.org/10.26686/wgtn.26421523
"recovery_scale"
