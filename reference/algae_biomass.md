# algae Biomass

The algae resource pool represents a combination of algae turf mats,
macroalgae, and the epilithic algae matrix (not including detritus). It
is not size structured to reflect the fact that herbivorous fish feed on
algae regardless of their body size.

## Usage

``` r
algae_biomass(params)
```

## Arguments

- params:

  MizerParams

## Value

The algae biomass in grams

## Details

A low standing biomass relative to production is ecologically expected
on a grazed reef: algal turfs are highly productive but are cropped down
by herbivores, keeping standing biomass low even though production is
high (see
[`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md)).
This is well documented empirically – for example Carpenter, R.C.
(1990). Mass mortality of *Diadema antillarum*: I. Long-term effects on
sea urchin population dynamics and coral reef algal communities. Marine
Biology, 104, 67-77, which found that algal cover roughly doubled within
months once a dominant grazer's population collapsed, with grazing
pressure otherwise removed. Note that this ecological expectation is
about the low value being *relative to a realistic absolute scale* (the
literature range is roughly 10-100 g dry weight/m^2, see
[`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)'s
`algae_growth_initial` docs) – it does not by itself explain the much
smaller absolute values (~1e-9 g/m^2) currently reported for the bundled
models, which is a separate, known calibration-scale issue (see
`inst/to-do-list.txt`'s "FUTURE WORK PLAN" entry).
