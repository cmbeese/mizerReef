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
#' For details on parameter derivation, see Chapter 3 of the PhD Thesis:
#' Beese, C. (2025). "Modelling Coral Reef Futures: Exploring the role of
#' structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington.
#'
#' Recalibrated in 2026 to follow the structure of Rogers et al. (2018). The
#' table also holds the calibrated values that `caribbean_3_model` is built
#' from: `w_min` (fish recruit at 10^-1.5 g), `gamma`, `h`, `ks`, `q`,
#' `erepro` and `R_max`, and the interactions with plankton, algae and
#' detritus. Only predators eat plankton; herbivores eat only algae and
#' invertebrates only detritus. Herbivores' `interaction_algae` and
#' invertebrates' `interaction_detritus` give an algae standing biomass of
#' 350 g/m^2 and a detritus lifetime of 3.85 hours. See
#' `inst/scripts/caribbean-3-calibration-clean.r`.
#'
#' @keywords datasets
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' http://theses.ncl.ac.uk/jspui/handle/10443/3574
#' @references Rogers, A., Blanchard, J.L., Newman, S.P., Dryden, C.S. &
#' Mumby, P.J. (2018). High refuge availability on coral reefs increases the
#' vulnerability of reef-associated predators to overexploitation. Ecology,
#' 99, 450-463. https://doi.org/10.1002/ecy.2103
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
#' groups: predators, herbivores, and invertebrates. Only predators eat other
#' groups. Their interactions were calibrated in 2026 so that adult
#' predators' diet is 33% predators, 17% herbivores and 50% invertebrates,
#' the proportions used by Rogers et al. (2018). Herbivores eat only algae and
#' invertebrates only detritus; those interactions are in
#' `caribbean_3_species`.
#'
#' @keywords datasets
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' @references Rogers, A., Blanchard, J.L., Newman, S.P., Dryden, C.S. &
#' Mumby, P.J. (2018). High refuge availability on coral reefs increases the
#' vulnerability of reef-associated predators to overexploitation. Ecology,
#' 99, 450-463. https://doi.org/10.1002/ecy.2103
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
#' herbivores, and invertebrates. It is built from `caribbean_3_species`,
#' `caribbean_3_interaction` and `karpata_refuge` by
#' `inst/scripts/caribbean-3-calibration-clean.r`, which projects the model
#' from its initial abundances to steady state. This object can be used
#' directly with mizer functions for simulation and analysis.
#'
#' Recalibrated in 2026 to follow the structure of Rogers et al. (2018):
#' herbivores eat only algae, invertebrates only detritus, and only predators
#' eat plankton (up to 10^-1.5 g). Fish compete for refuge holes from 2.5 cm
#' (`method = "competitive"`, hole densities from Karpata reef). Holes
#' outnumber fish, so about 98% of fish between 2.5 and 20 cm are protected.
#' Algae has a carrying capacity; detritus does not. At steady state:
#' * biomass above 3.125 g is 107 (predators), 34.1 (herbivores) and 40.1
#'   (invertebrates) g/m^2, matching the FORCE-survey targets (107/34/40);
#' * there are about 10 herbivore juveniles (0.03-1.6 g) per m^2;
#' * adult predators' diet is 33% predators, 17% herbivores and 50%
#'   invertebrates;
#' * predators mature at age 2 and herbivores at 1.7;
#' * algae standing biomass is 350 g/m^2 and the detritus lifetime 3.85 hours.
#'
#' For details on parameter derivation and implementation, see Chapter 3 of
#' the PhD Thesis:
#' Beese, C. (2025). "Modelling Coral Reef Futures: Exploring the role of
#' structural complexity in sustaining ecosystem services."
#' Victoria University of Wellington.
#'
#' @keywords datasets
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish community
#' structure. PhD Thesis, Newcastle University.
#' @references Rogers, A., Blanchard, J.L., Newman, S.P., Dryden, C.S. &
#' Mumby, P.J. (2018). High refuge availability on coral reefs increases the
#' vulnerability of reef-associated predators to overexploitation. Ecology,
#' 99, 450-463. https://doi.org/10.1002/ecy.2103
#' @format A [mizerReef-class] object.
#'
#' The object was created with [newReefParams()] and is stored in the package's
#' `data/` directory. R's standard lazy-loading preserves its S3 class vector
#' and extension metadata.
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
#' \code{biomass_observed} is NA for \code{eels}, \code{pred_crypt},
#' \code{farm_damsel}, and \code{inverts} -- the thesis's own visual-survey
#' methodology excluded fish under 10cm and struggled to capture
#' small-bodied/cryptic groups, so no FORCE observation exists for these four
#' (Beese 2025, Chapter 4). Several parameters were revised in 2026 during
#' recalibration against a corrected senescence-mortality formula (see
#' \code{\link{caribbean_10_model}}'s documentation and
#' \code{inst/scripts/Caribbean_10_model-calibration.R}'s design note for the
#' full reasoning and citations): \code{satiation} for \code{parrotfish},
#' \code{farm_damsel}, and \code{herbs} (FALSE -> TRUE, activating literature/
#' thesis-derived intake caps that were previously gated off); \code{age_mat}
#' for \code{parrotfish} (2 -> 1.6, Rivera Hernandez & Shervette 2025) and
#' \code{pred_grab} (2 -> 5.4, derived from the thesis's own stored growth
#' curve combined with an independently-reported maturity length for
#' \emph{Lutjanus apodus}, FishBase); and \code{biomass_observed} for
#' \code{eels} and \code{farm_damsel} (NA -> soft literature-derived targets
#' of 10 and 0.4 g/m^2 respectively, Gilbert et al. 2005; Vermeij et al.
#' 2015 -- these are not FORCE observations and are not treated as exact).
#'
#' @keywords datasets
#' @seealso \code{\link{caribbean_10_model}}, \code{\link{caribbean_10_interaction}},
#'   \code{\link{karpata_refuge}}
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' @references Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
#'   https://doi.org/10.26686/wgtn.26421523
#' @references Rivera Hernandez, J.M. & Shervette, V.R. (2025). Addressing
#' life history information gaps for Caribbean parrotfishes: queen parrotfish
#' Scarus vetula and stoplight parrotfish Sparisoma viride. Environmental
#' Biology of Fishes, 108, 179-198. https://doi.org/10.1007/s10641-024-01651-x
#' @references Gilbert, M., Rasmussen, J.B. & Kramer, D.L. (2005). Estimating
#' the density and biomass of moray eels (Muraenidae) using a modified visual
#' census method for hole-dwelling reef fauna. Environmental Biology of
#' Fishes, 73, 415-426. https://doi.org/10.1007/s10641-005-2228-2
#' @references Vermeij, M.J.A., DeBey, H., Grimsditch, G., Brown, J., Obura,
#' D., DeLeon, R. & Sandin, S.A. (2015). Negative effects of gardening
#' damselfish Stegastes planifrons on coral health depend on predator
#' abundance. Marine Ecology Progress Series, 528, 289-296.
#' https://doi.org/10.3354/meps11243
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
#' calibrated to match observed biomasses and built with mizerReef 2.0.0 /
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
#' Recalibrated in 2026 against a corrected senescence-mortality formula
#' (\code{reefSenMort()}), the same underlying issue that required
#' recalibrating \code{\link{caribbean_3_model}}. Biomass, growth, and
#' reproduction were all re-tuned (verified: repeated \code{reefSteady()}
#' calls leave biomass unchanged). Final biomass is within 79-121\% of the
#' six FORCE-observed targets (\code{pred_eng}, \code{pred_grab},
#' \code{pred_inv}, \code{pred_plank}, \code{parrotfish}, \code{herbs});
#' growth (age at maturity) matches its target closely for five of those six,
#' and reasonably (within a factor of ~3) for \code{pred_grab} after
#' correcting its own \code{age_mat} target (see
#' \code{\link{caribbean_10_species}}'s documentation). \code{eels} and
#' \code{farm_damsel} -- which have no FORCE-observed biomass target -- were
#' given soft, literature-derived targets and brought up from near-zero to a
#' modest, non-vanishing presence, deliberately not matched exactly (doing so
#' costs disproportionate accuracy on the six well-observed groups). A
#' version of this model that ignores \code{eels}/\code{farm_damsel} entirely
#' and matches the six FORCE-observed groups more tightly is archived (not
#' bundled as a package dataset, to avoid doubling the recalibration burden)
#' at \code{inst/archive/caribbean_10_model_tight_biomass_fit.rds} -- load
#' with \code{mizer::readParams()}. The complexity-increases-biomass
#' relationship central to the associated manuscript/thesis chapter was
#' explicitly re-verified against this recalibrated model (competitive vs.
#' non-complex refuge), not merely assumed to still hold. See
#' \code{inst/scripts/Caribbean_10_model-calibration.R}'s design note for the
#' full reasoning, including why growth-matching required tuning one species
#' at a time rather than mizerExperimental's usual batched workflow.
#'
#' @keywords datasets
#' @seealso \code{\link{caribbean_10_species}}, \code{\link{caribbean_10_interaction}},
#'   \code{\link{karpata_refuge}}, \code{\link{newReefParams}}
#' @references FORCE dataset. From: Dryden, C. (2017). Habitat structural
#' complexity of Caribbean coral reefs and its relationships with fish
#' community structure. PhD Thesis, Newcastle University.
#' @format A [mizerReef-class] object with 10 species, competitive refuge, and
#'   algae/detritus components.
#'
#' The object was created with [newReefParams()] and is stored in the package's
#' `data/` directory. R's standard lazy-loading preserves its S3 class vector
#' and extension metadata.
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
