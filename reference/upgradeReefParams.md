# Upgrade a mizerReef 1.x params object to the current (2.0.0+) layout

`mizerReef` 2.x replaced the direct rate-function overrides and flat
`other_params` layout used by mizerReef 1.x (the version accompanying
the original thesis, on this package's `main` branch) with S3 dispatch
on a `mizerReef` marker class and a nested `other_params` layout.
Objects created with mizerReef 1.x are never automatically upgraded by
[`mizer::validParams()`](https://sizespectrum.org/mizer/reference/validParams.html),
because they were never registered as a `"mizerReef"` extension in the
first place – there is nothing for the automatic
[`upgrade.mizerReef()`](https://cmbeese.github.io/mizerReef/reference/upgrade.mizerReef.md)
hook to trigger on. Call this function once, by hand, on any params
object created with mizerReef 1.x to bring it up to the current
structural layout.

## Usage

``` r
upgradeReefParams(params)
```

## Arguments

- params:

  A `MizerParams` object created with mizerReef 1.x (i.e. from code on
  this package's `main`/thesis branch, not a `mizerReef` object already
  created with this version of the package).

## Value

A `mizerReef` object with the current structural layout.

## Details

This is a **structural** migration only. Tuned values that came from
your original model – `algae_growth`, `detritus_external`, capacities,
the refuge method and its parameters, species parameters – are carried
over as they were, not reset to package defaults and not silently
re-tuned. In particular, the refuge profile
(`refuge`/`bin.id`/`refuge_lengths`) is recomputed from your preserved
method/method_params/`a_bar`/`b_bar` using the current
[`getRefuge()`](https://cmbeese.github.io/mizerReef/reference/getRefuge.md)
logic, since that is a deterministic function of already-preserved
inputs, not a recalibration.

Because your original algae/detritus tuning is carried over unchanged,
an upgraded object is not guaranteed to be at a genuine joint steady
state under the current algae/detritus mechanism (see
[`vignette("extension_mechanism")`](https://cmbeese.github.io/mizerReef/articles/extension_mechanism.md)
and the algae/detritus sections of
[`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)/[`setDetritusParams()`](https://cmbeese.github.io/mizerReef/reference/setDetritusParams.md)
for what changed). Run
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
on the result if you want a fresh steady state under the current
mechanism.

## See also

[`upgrade.mizerReef()`](https://cmbeese.github.io/mizerReef/reference/upgrade.mizerReef.md)
for the automatic, incremental micro-version upgrades that run every
time
[`mizer::validParams()`](https://sizespectrum.org/mizer/reference/validParams.html)
is called on an already-registered `mizerReef` object – this function is
for the one-off 1.x -\> 2.x jump only.

## Examples

``` r
if (FALSE) { # \dontrun{
# `old_params` is a MizerParams object created with mizerReef 1.x (the
# version on this package's `main`/thesis branch), e.g. loaded from a
# saved .rds/.rda file from before upgrading the package.
params <- upgradeReefParams(old_params)
} # }
```
