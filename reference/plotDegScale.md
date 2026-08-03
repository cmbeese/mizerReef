# Compare built-in degradation trajectories as faceted heatmaps

Creates a faceted heatmap comparing the three built-in degradation
trajectories (rubble, algae, and recovery) side by side. Each panel
shows how refuge density scaling factors change over time and across
refuge size bins following a bleaching event.

## Usage

``` r
plotDegScale(return_data = FALSE)
```

## Arguments

- return_data:

  Logical. If TRUE, returns the combined long-format data frame instead
  of the plot. Default FALSE.

## Value

A ggplot2 object, or a data frame if `return_data = TRUE`.

## Details

The colour scale shows scaling factors relative to pre-disturbance
refuge density: values below 1 (red/pink) indicate a reduction, 1 (grey)
indicates no change, and values above 1 (green/blue) indicate an
increase.

## See also

[`plotDegradationScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegradationScale.md),
[`setDegradation()`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md),
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md),
[rubble_scale](https://cmbeese.github.io/mizerReef/reference/rubble_scale.md),
[algae_scale](https://cmbeese.github.io/mizerReef/reference/algae_scale.md),
[recovery_scale](https://cmbeese.github.io/mizerReef/reference/recovery_scale.md)

Other plotting functions:
[`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
[`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md),
[`plotDegradationScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegradationScale.md),
[`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
[`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md),
[`plotRefugeDensity()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeDensity.md),
[`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md),
[`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md),
[`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md),
[`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md),
[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md),
[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md),
[`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
