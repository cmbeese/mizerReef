# Competitive method refuge parameters for Aquarius reef research station

This data frame contains refuge density parameters for the competitive
refuge method. It records the density of refuges (no./m^2) in each of
ten 5 cm fish-length bins (0-50 cm), derived from field measurements at
the Aquarius reef research station (FORCE dataset, Dryden 2017), from
the same field hole-density survey pipeline as
[`karpata_refuge`](https://cmbeese.github.io/mizerReef/reference/karpata_refuge.md).
Like `karpata_refuge`, this profile is general-purpose and can be used
with any mizerReef model, not just one particular example model.

## Usage

``` r
aquarius_refuge
```

## Format

A data frame with 10 rows and 3 columns: `start_L` (start of length bin
in cm), `end_L` (end of length bin in cm), and `refuge_density` (refuge
density in no./m^2).

## Source

FORCE dataset, Aquarius reef research station (Dryden 2017)

## Details

Use this with `method = "competitive"` in
[`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md).

## References

FORCE dataset. From: Dryden, C. (2017). Habitat structural complexity of
Caribbean coral reefs and its relationships with fish community
structure. PhD Thesis, Newcastle University.

Beese, C. (2025). PhD Thesis. "Modelling Coral Reef Futures: Exploring
the role of structural complexity in sustaining ecosystem services."
Victoria University of Wellington.
https://doi.org/10.26686/wgtn.26421523

## See also

[`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md),
[`karpata_refuge`](https://cmbeese.github.io/mizerReef/reference/karpata_refuge.md)
