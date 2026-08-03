# species_params dataframe for a simple, generic Caribbean reef with 3 species groups

This data frame contains species-level parameters for three functional
groups: predators, herbivores, and invertebrates. Parameters include
biomasses from empirical data collected in Bonaire (FORCE dataset, see:
Dryden, C. (2017), "Habitat structural complexity of Caribbean coral
reefs and its relationships with fish community structure", Newcastle
University), as well as growth rates, maturation ages, and other
life-history traits derived from FishBase (www.fishbase.org) and
published literature. See thesis for further details on parameter
sources and derivation.

## Usage

``` r
caribbean_3_species
```

## Format

data frame

## Source

PhD Thesis, FORCE dataset

## Details

For details on parameter derivation and model calibration, see Chapter 3
of the PhD Thesis: Beese, C. (2025). "Modelling Coral Reef Futures:
Exploring the role of structural complexity in sustaining ecosystem
services." Victoria University of Wellington.

Two herbivore parameters were revised in 2026 during recalibration
against a corrected senescence-mortality formula (see
`caribbean_3_model`'s documentation and
`inst/scripts/Caribbean_3_model-calibration.R`'s design note):
`satiation` (FALSE -\> TRUE) and `age_mat` (4 -\> 1.6, the age at median
sexual maturity reported for the stoplight parrotfish, *Sparisoma
viride*, by Rivera Hernandez & Shervette (2025), replacing a value that
appears to have conflated it with that species' age at sexual
transition).

## References

FORCE dataset. From: Dryden, C. (2017). Habitat structural complexity of
Caribbean coral reefs and its relationships with fish community
structure. PhD Thesis, Newcastle University.
http://theses.ncl.ac.uk/jspui/handle/10443/3574

Rivera Hernandez, J.M. & Shervette, V.R. (2025). Addressing life history
information gaps for Caribbean parrotfishes: queen parrotfish Scarus
vetula and stoplight parrotfish Sparisoma viride. Environmental Biology
of Fishes, 108, 179-198. https://doi.org/10.1007/s10641-024-01651-x

Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures: Exploring
the role of structural complexity in sustaining ecosystem services."
Victoria University of Wellington.
https://doi.org/10.26686/wgtn.26421523
