# Species parameters for a generic Caribbean reef (10 functional groups)

A species parameter data frame for a 10-group size-spectrum model of a
Caribbean coral reef. Functional groups include multiple predator guilds
(engulf, grab, eel-like, cryptic, invertebrate, planktivorous),
parrotfish, farming damselfish, other herbivores, and invertebrates.
Parameters are derived from field measurements at Karpata Reef, Bonaire
(FORCE dataset, Dryden 2017) and published life-history data (FishBase).

## Usage

``` r
caribbean_10_species
```

## Format

A data frame with one row per functional group and columns for
life-history parameters (w_max, w_mat, age_mat, h, beta, sigma, etc.),
refuge use (`refuge_user`, `blocked_pred`, `satiation`), and resource
interactions (`interaction_algae`, `interaction_detritus`).

## Source

PhD Thesis, FORCE dataset

## Details

`biomass_observed` is NA for `eels`, `pred_crypt`, `farm_damsel`, and
`inverts` – the thesis's own visual-survey methodology excluded fish
under 10cm and struggled to capture small-bodied/cryptic groups, so no
FORCE observation exists for these four (Beese 2025, Chapter 4). Several
parameters were revised in 2026 during recalibration against a corrected
senescence-mortality formula (see
[`caribbean_10_model`](https://cmbeese.github.io/mizerReef/reference/caribbean_10_model.md)'s
documentation and `inst/scripts/Caribbean_10_model-calibration.R`'s
design note for the full reasoning and citations): `satiation` for
`parrotfish`, `farm_damsel`, and `herbs` (FALSE -\> TRUE, activating
literature/ thesis-derived intake caps that were previously gated off);
`age_mat` for `parrotfish` (2 -\> 1.6, Rivera Hernandez & Shervette
2025) and `pred_grab` (2 -\> 5.4, derived from the thesis's own stored
growth curve combined with an independently-reported maturity length for
*Lutjanus apodus*, FishBase); and `biomass_observed` for `eels` and
`farm_damsel` (NA -\> soft literature-derived targets of 10 and 0.4
g/m^2 respectively, Gilbert et al. 2005; Vermeij et al. 2015 – these are
not FORCE observations and are not treated as exact).

## References

FORCE dataset. From: Dryden, C. (2017). Habitat structural complexity of
Caribbean coral reefs and its relationships with fish community
structure. PhD Thesis, Newcastle University.

Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
https://doi.org/10.26686/wgtn.26421523

Rivera Hernandez, J.M. & Shervette, V.R. (2025). Addressing life history
information gaps for Caribbean parrotfishes: queen parrotfish Scarus
vetula and stoplight parrotfish Sparisoma viride. Environmental Biology
of Fishes, 108, 179-198. https://doi.org/10.1007/s10641-024-01651-x

Gilbert, M., Rasmussen, J.B. & Kramer, D.L. (2005). Estimating the
density and biomass of moray eels (Muraenidae) using a modified visual
census method for hole-dwelling reef fauna. Environmental Biology of
Fishes, 73, 415-426. https://doi.org/10.1007/s10641-005-2228-2

Vermeij, M.J.A., DeBey, H., Grimsditch, G., Brown, J., Obura, D.,
DeLeon, R. & Sandin, S.A. (2015). Negative effects of gardening
damselfish Stegastes planifrons on coral health depend on predator
abundance. Marine Ecology Progress Series, 528, 289-296.
https://doi.org/10.3354/meps11243

## See also

[`caribbean_10_model`](https://cmbeese.github.io/mizerReef/reference/caribbean_10_model.md),
[`caribbean_10_interaction`](https://cmbeese.github.io/mizerReef/reference/caribbean_10_interaction.md),
[`karpata_refuge`](https://cmbeese.github.io/mizerReef/reference/karpata_refuge.md)
