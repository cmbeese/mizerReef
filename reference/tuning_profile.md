# Constant refuge profile for tuning steady states

This data frame provides a constant refuge profile for model calibration
and steady-state tuning. It contains start and end lengths for size bins
and sets `prop_protect` to 60% for all bins up to 50 cm in length.

## Usage

``` r
tuning_profile
```

## Format

data frame

## Source

Beese PhD Thesis

## Details

These parameters are intended to be used with the "binned" method, which
is independent of density. The tuning profile should be used while
calibrating biomass and growth rates before switching to the competitive
method. See the "Getting Started" vignette for more details.

## References

Beese, C. (2025). PhD Thesis. Victoria University of Wellington.
<https://doi.org/10.26686/thesis.123456>
