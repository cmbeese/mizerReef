# Scale background down by a factor

**\[superseded\]** Use
[`scaleDownPlankton()`](https://cmbeese.github.io/mizerReef/reference/scaleDownPlankton.md)
instead, which leaves algae and detritus completely unchanged and keeps
each species' reproduction level.

## Usage

``` r
scaleReefBackground(params, factor)
```

## Arguments

- params:

  a mizer model object

- factor:

  A number giving the factor by which the background abundance will be
  reduced

## Details

Multiplies the fish abundances by `factor` with
[`scaleReefAbundance()`](https://cmbeese.github.io/mizerReef/reference/scaleReefAbundance.md),
which resets every species' reproduction level to 1/2, and then scales
the whole model by `1 / factor` with
[`mizer::scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.html).
Algae and detritus biomasses end up divided by `factor`, and their
encounter coefficients `rho` multiplied by `factor`.
