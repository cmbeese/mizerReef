# Scale down the plankton resource by a factor

Reduces the abundance of the size-structured plankton resource by
`factor` relative to the fish, algae and detritus. Every species then
encounters `factor` times as much fish prey relative to plankton, which
shifts predators' diets away from plankton and towards fish during
calibration.

## Usage

``` r
scaleDownPlankton(params, factor)
```

## Arguments

- params:

  A MizerParams object

- factor:

  A number greater than 0 giving the factor by which the plankton
  abundance is reduced

## Value

An updated MizerParams object

## Details

The plankton abundance, its carrying capacity and `kappa` are divided by
`factor`, and every species' search volume and `gamma` are multiplied by
`factor`. So at first each species encounters as much plankton as
before, and `factor` times as much fish prey. Algae and detritus are
left completely unchanged, so herbivores and invertebrates still
encounter as much of them as before.

The fish then graze the plankton `factor` times harder, though. Once the
model is run back to steady state, the plankton settles below 1/`factor`
of its former abundance wherever grazing is strong compared with the
plankton's regrowth rate, so the eventual shift towards fish prey is
larger than `factor` alone suggests.

Background species (`is_background = TRUE`) count as background food, as
in
[`mizerExperimental::scaleDownBackground()`](https://sizespectrum.org/mizerExperimental/reference/scaleDownBackground.html),
and are scaled down with the plankton. Size-structured resources added
by other extensions are not scaled.

If the model uses Beverton-Holt recruitment, each species keeps its
reproduction level: `erepro` and `R_max` are retuned with
[`mizer::setBevertonHolt()`](https://sizespectrum.org/mizer/reference/setBevertonHolt.html).
Any other recruitment function is left as it is.

This changes the model's rates, so the model is no longer at steady
state. Follow it with
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
and then re-match biomasses and growth, for example with
[`mizer::matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.html)
and
[`matchReefGrowth()`](https://cmbeese.github.io/mizerReef/reference/matchReefGrowth.md).

[`mizerExperimental::scaleDownBackground()`](https://sizespectrum.org/mizerExperimental/reference/scaleDownBackground.html)
also scales down the plankton, but it resets every species' reproduction
level to 1/4, and it divides the algae and detritus biomasses by
`factor`. It multiplies their encounter coefficients `rho` by `factor`
to make up for that, which leaves the encounter rates unchanged but
makes algae and detritus turn over `factor` times faster. It also
divides the external detritus flux by `factor`, so the detritus budget
no longer balances until
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
retunes it.

## See also

[`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md)

## Examples

``` r
data(caribbean_3_model)
params <- scaleDownPlankton(caribbean_3_model, factor = 2)
```
