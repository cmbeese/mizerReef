# Degradation scaling matrix for the rubble trajectory

A 10 x 15 numeric matrix of multiplicative scaling factors for refuge
density under the rubble disturbance trajectory. Each row corresponds to
a 5 cm refuge size bin (0-5 cm through 45-50 cm) and each column to a
simulation year (column 1 = bleaching year, columns 2-15 =
post-bleaching years 1-14). Values are applied as
`new_rd = scale * old_rd` inside
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md),
so a value of 1.4 means a 40\\ refuge density and 0.4 means a 60\\

## Usage

``` r
rubble_scale
```

## Format

A numeric matrix with 10 rows (refuge size bins) and 15 columns (time
steps). Row names are size bin labels (e.g. `"0 to 5"`); column names
are integers 1 to 15.

## Details

In the rubble trajectory, large structural refuges collapse while an
initial pulse of rubble temporarily increases the availability of small
refuges (0-5 cm bin) before all sizes return to baseline.

Values are derived from `inst/data-csv/deg_scales.csv` by computing
`multiplier = 1 + delta` for each cell.

## References

Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
https://doi.org/10.26686/wgtn.26421523

## See also

[`algae_scale`](https://cmbeese.github.io/mizerReef/reference/algae_scale.md),
[`recovery_scale`](https://cmbeese.github.io/mizerReef/reference/recovery_scale.md),
[`setDegradation`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md),
[`reefDegrade`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
