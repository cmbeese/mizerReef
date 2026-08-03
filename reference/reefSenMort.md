# Expanding external mortality rate to include senescence

Senescence mortality represents mortality caused by background sources
such as illness or age, in addition to residual natural mortality (see
[`setExtMortParams()`](https://cmbeese.github.io/mizerReef/reference/setExtMortParams.md)'s
"Residual natural mortality" section), which represents all other
mortality that is not due to fishing or predation by predators included
in the model.

## Usage

``` r
reefSenMort(params, ...)
```

## Arguments

- params:

  A MizerParams object

- ...:

  Unused

## Value

A named two dimensional array (species x size) with the senescence
mortality rates.

## Details

Users can change `sen_prop`/`sen_curve` via
[`setExtMortParams()`](https://cmbeese.github.io/mizerReef/reference/setExtMortParams.md).

## Senescence mortality

     Senescence mortality \eqn{\mu_{sen.i}(w)} is used to represent
     mortality caused by background sources such as illness or age. The
     rate of senescence mortality (in 1/year) is given by:

     \deqn{\mu_{sen.i}(w) = \mathtt{sen\_prop}\,
                                 \left[\max\left(0,\;
                                 \frac{\log_{10}(w)}{\log_{10}(w_{max.i})}
                                 \right)\right]^{\mathtt{sen\_curve}}}{
              \mu_{sen.i}(w) = sen_prop *
              max(0, log10(w)/log10(w_{max.i}))^sen_curve}

     where \eqn{\mathtt{sen\_curve}} is the exponent shaping the
     senescence curve and \eqn{\mathtt{sen\_prop}} is the rate the curve
     approaches as \eqn{w \to w_{max.i}} (where the ratio is exactly 1).
     The ratio is floored at zero before being raised to the
     \eqn{\mathtt{sen\_curve}} power, since it is negative for individuals
     below 1 gram (where \eqn{\log_{10}(w) < 0}), which would otherwise
     raise a negative number to a fractional power -- those individuals
     get exactly zero senescence mortality.
