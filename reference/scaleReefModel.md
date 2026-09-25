# Scale model parameters

This function scales various model parameters by a given factor.

## Usage

``` r
scaleReefModel(params, factor)

# S3 method for class 'mizerReef'
scaleModel(params, factor, ...)
```

## Arguments

- params:

  a mizer model object

- factor:

  a numeric value by which to scale the model

- ...:

  Unused

## Value

a mizer model object with scaled parameters

## Details

The fish abundances, the plankton and the algae and detritus biomasses
are multiplied by `factor`, and the search volume and the algae and
detritus encounter coefficients `rho` are divided by it, so every
encounter rate is unchanged at first. The external detritus flux is
multiplied by `factor`.

Two things are properties of the reef, not of the fish, and are
deliberately not scaled. The algae production rate is a literature
value, so algae biomass then settles to whatever level the new grazing
pressure allows. The refuge density of the competitive refuge method
comes from data, so scaling the fish abundance changes how many fish
find a refuge. The algae and detritus carrying capacities (with
`use_UR_cc = TRUE`) are not scaled either.

On a mizerReef model, mizer's own
[`mizer::scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.html)
does exactly the same as `scaleReefModel()`. So mizer and
mizerExperimental functions that rescale a model through `scaleModel()`,
such as
[`mizer::calibrateBiomass()`](https://sizespectrum.org/mizer/reference/calibrateBiomass.html)
and
[`mizerExperimental::scaleDownBackground()`](https://sizespectrum.org/mizerExperimental/reference/scaleDownBackground.html),
scale algae and detritus in the same way.

## See also

[`scaleDownPlankton()`](https://cmbeese.github.io/mizerReef/reference/scaleDownPlankton.md)
