# Expanding external mortality rate to include senescence

Expanding external mortality rate to include senescence

## Usage

``` r
reefSenMort(params, n, n_pp, n_other, t, ...)
```

## Arguments

- params:

  A MizerParams object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components of the ecosystem

- t:

  The time for which to do the calculation (Not used by standard mizer
  rate functions but useful for extensions with time-dependent
  parameters.)

- ...:

  Unused

## Value

A named two dimensional array (species x size) with the senescence
mortality rates.

## Senescence mortality

     Senescence mortality \eqn{\mu_{sen.i}(w)} is used to represent
     mortality caused by external sources such as illness or age. This is
     addition to external mortality, \eqn{\mu_{ext.i}(w)}, which represents
     all mortality that is not due to fishing or predation by predators
     included in the model. The rate of senescence mortality is given by:

     \deqn{\mu_{sen.i}(w) = k_{sen}\left(\frac{w}{w_{max.i}}\right)^{p_{sen}}}
          {\mu_{sen.i}(w) = k_{sen} (w/w_{max.i})^{p_{sen}}

     where \eqn{k_{sen}} and \eqn{p_sen} are constants defining the shape
     of the senescence curve and \eqn{w_max.i} is the maximum size
     of species \eqn{i} in grams.

     Users can change all constants with the `setSenMortParams()` function.
