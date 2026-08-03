# Get the size specific senescence mortality rate

Returns the rate of senescence mortality at each size by functional
group.

## Usage

``` r
getSenMort(
  params,
  n = initialN(params),
  n_pp = params@initial_n_pp,
  n_other = initialNOther(params),
  t = 0,
  ...
)
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

## See also

Other rate functions:
[`getDegrade()`](https://cmbeese.github.io/mizerReef/reference/getDegrade.md),
[`getVulnerable()`](https://cmbeese.github.io/mizerReef/reference/getVulnerable.md)
