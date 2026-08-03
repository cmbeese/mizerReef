# mizerReef model for a generic Caribbean reef (10 functional groups)

A calibrated `mizerReef` object for a generic Caribbean coral reef with
ten functional groups, parameterised from field data collected at
Karpata Reef, Bonaire (FORCE dataset, Dryden 2017). Functional groups
span six predator guilds, parrotfish, farming damselfish, other
herbivores, and invertebrates. The model includes competitive predation
refuge (using
[`karpata_refuge`](https://cmbeese.github.io/mizerReef/reference/karpata_refuge.md))
and benthic resources (algae and detritus).

## Usage

``` r
caribbean_10_model
```

## Format

A `mizerReef` object (extends
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html))
with 10 species, competitive refuge, and algae/detritus components.

## Source

PhD Thesis, FORCE dataset

## Details

Built from
[`caribbean_10_species`](https://cmbeese.github.io/mizerReef/reference/caribbean_10_species.md),
[`caribbean_10_interaction`](https://cmbeese.github.io/mizerReef/reference/caribbean_10_interaction.md),
and
[`karpata_refuge`](https://cmbeese.github.io/mizerReef/reference/karpata_refuge.md);
calibrated to match observed biomasses and built with mizerReef 2.0.0 /
mizer 3.1.0. Detritus resource dynamics are tuned via
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md) so
that they are genuinely at steady state (`dB/dt = 0`) for the model's
calibrated abundances. Algae biomass is likewise tuned by
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md) to
its own steady state, but for a fixed, literature-informed production
rate (`algae_growth`, left at
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)'s
default of `2000` grams per square meter per year – see
[`setAlgaeParams`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md))
that is not itself retuned to match consumption, since real algal
production is not driven by grazer demand; an earlier untuned version of
this object (algae biomass not yet solved for its steady state) is
archived at `inst/archive/caribbean_10_model_untuned.rda` for reference.

Recalibrated in 2026 against a corrected senescence-mortality formula
([`reefSenMort()`](https://cmbeese.github.io/mizerReef/reference/reefSenMort.md)),
the same underlying issue that required recalibrating
[`caribbean_3_model`](https://cmbeese.github.io/mizerReef/reference/caribbean_3_model.md).
Biomass, growth, and reproduction were all re-tuned (verified: repeated
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
calls leave biomass unchanged). Final biomass is within 79-121\\ six
FORCE-observed targets (`pred_eng`, `pred_grab`, `pred_inv`,
`pred_plank`, `parrotfish`, `herbs`); growth (age at maturity) matches
its target closely for five of those six, and reasonably (within a
factor of ~3) for `pred_grab` after correcting its own `age_mat` target
(see
[`caribbean_10_species`](https://cmbeese.github.io/mizerReef/reference/caribbean_10_species.md)'s
documentation). `eels` and `farm_damsel` – which have no FORCE-observed
biomass target – were given soft, literature-derived targets and brought
up from near-zero to a modest, non-vanishing presence, deliberately not
matched exactly (doing so costs disproportionate accuracy on the six
well-observed groups). A version of this model that ignores
`eels`/`farm_damsel` entirely and matches the six FORCE-observed groups
more tightly is archived (not bundled as a package dataset, to avoid
doubling the recalibration burden) at
`inst/archive/caribbean_10_model_tight_biomass_fit.rds` – load with
[`mizer::readParams()`](https://sizespectrum.org/mizer/reference/saveParams.html).
The complexity-increases-biomass relationship central to the associated
manuscript/thesis chapter was explicitly re-verified against this
recalibrated model (competitive vs. non-complex refuge), not merely
assumed to still hold. See
`inst/scripts/Caribbean_10_model-calibration.R`'s design note for the
full reasoning, including why growth-matching required tuning one
species at a time rather than mizerExperimental's usual batched
workflow.

## References

FORCE dataset. From: Dryden, C. (2017). Habitat structural complexity of
Caribbean coral reefs and its relationships with fish community
structure. PhD Thesis, Newcastle University.

Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
https://doi.org/10.26686/wgtn.26421523

## See also

[`caribbean_10_species`](https://cmbeese.github.io/mizerReef/reference/caribbean_10_species.md),
[`caribbean_10_interaction`](https://cmbeese.github.io/mizerReef/reference/caribbean_10_interaction.md),
[`karpata_refuge`](https://cmbeese.github.io/mizerReef/reference/karpata_refuge.md),
[`newReefParams`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
