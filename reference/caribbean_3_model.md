# mizerParams object for a simple, generic Caribbean reef with 3 species groups

A `MizerParams` object for a simple Caribbean reef ecosystem model with
three functional groups: predators, herbivores, and invertebrates. It is
built using the `caribbean_3_species` and `caribbean_3_interaction`
data, and includes all parameter settings and calibrations as described
in the associated PhD thesis and FORCE dataset. This object can be used
directly with mizer functions for simulation and analysis.

## Usage

``` r
caribbean_3_model
```

## Format

mizerParams object

## Source

PhD Thesis, FORCE dataset

## Details

Recalibrated in 2026 against a corrected senescence-mortality formula
([`reefSenMort()`](https://cmbeese.github.io/mizerReef/reference/reefSenMort.md))
– fish abundances, growth/consumption rates and reproduction parameters
were all re-tuned to restore a genuine steady state under the corrected
code (verified: repeated
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
calls leave biomass unchanged). Two herbivore species_params were also
changed from their original thesis values as part of this: `satiation`
is now `TRUE` (was `FALSE`) since herbivore biomass has no
density-dependent brake without some cap on individual intake once
mortality is realistically low – the realised feeding level is
consistently close to 1, still consistent with the citations behind the
package's usual herbivore default (Caribbean herbivores' guts observed
to be full nearly continuously); and `age_mat` is now `1.6` (was `4`),
the age at median sexual maturity for *Sparisoma viride* (Rivera
Hernandez & Shervette 2025), replacing a value that appears to have
conflated it with that species' age at sexual transition. See
`inst/scripts/Caribbean_3_model-calibration.R`'s design note and
`vignettes/caribbean_3_model-description.Rmd` for the full reasoning.
Biomasses (predators/herbivores/inverts, g/m^2) are now 107.9/34.1/40.1,
matching the FORCE-survey targets (107/34/40) essentially exactly.

For details on parameter derivation, model calibration, and
implementation, see Chapter 3 of the PhD Thesis: Beese, C. (2025).
"Modelling Coral Reef Futures: Exploring the role of structural
complexity in sustaining ecosystem services." Victoria University of
Wellington.

## References

FORCE dataset. From: Dryden, C. (2017). Habitat structural complexity of
Caribbean coral reefs and its relationships with fish community
structure. PhD Thesis, Newcastle University.

Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures: Exploring
the role of structural complexity in sustaining ecosystem services."
Victoria University of Wellington.
https://doi.org/10.26686/wgtn.26421523
