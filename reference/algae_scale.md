# Degradation scaling matrix for the algae trajectory

A 10 x 15 numeric matrix of multiplicative scaling factors for refuge
density under the algae disturbance trajectory. Each row corresponds to
a 5 cm refuge size bin (0-5 cm through 45-50 cm) and each column to a
simulation year (column 1 = bleaching year, columns 2-15 =
post-bleaching years 1-14).

## Usage

``` r
algae_scale
```

## Format

A numeric matrix with 10 rows (refuge size bins) and 15 columns (time
steps). Row names are size bin labels; column names are integers 1 to
15.

## Details

In the algae trajectory, moderate coral loss is partially offset by
algal colonisation of rubble. Small refuges (0-5 cm bin) show a slight
overshoot above baseline (1.02x) during years 6-9 before returning to
pre-disturbance levels.

Values are derived from `inst/data-csv/deg_scales.csv` by computing
`multiplier = 1 + delta` for each cell.

## References

Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
https://doi.org/10.26686/wgtn.26421523

## See also

[`rubble_scale`](https://cmbeese.github.io/mizerReef/reference/rubble_scale.md),
[`recovery_scale`](https://cmbeese.github.io/mizerReef/reference/recovery_scale.md),
[`setDegradation`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md),
[`reefDegrade`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
