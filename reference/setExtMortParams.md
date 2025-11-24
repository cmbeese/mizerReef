# Set the parameters for external mortality

mizerReef models contain two sources of external mortality, residual
natural mortality and senescence. This function checks if given
parameters are valid, sets defaults, and then stores them in the
mizerParams object.

## Usage

``` r
setExtMortParams(params, ext_mort_params = NULL)
```

## Arguments

- params:

  MizerParams object

  Optional parameters:

- ext_mort_params:

  Named list containing desired mortality parameters

## Value

MizerParams object with updated mortality parameters

## Details

External mortality is implemented by default in mizerReef. You do not
need to set these parameters. This function should only be used to
change default values.

## Residual natural mortality

     Residual natural mortality accounts for any external predation or
     fishing mortality that is not explicitly included in the model. It is
     assumed to decrease allometrically with body size. Residual natural
     mortality is a rate with units 1/year given by:

     \deqn{\mu_{nat.i}(w) = \mu_{nat}\, w^{1-n}.}
          {\mu_{nat.i}(w) = \mu_{nat}\, w^{1-n}.}

      Here \eqn{\mu_{nat}} is the residual natural mortality rate at size
      1 g and \eqn{n} is the allometric scaling exponent. In mizerReef,
      these default to \eqn{\mu_{nat} = 0.2} and \eqn{n = 0.75}.

## Senescence mortality

     Senescence mortality \eqn{\mu_{sen.i}(w)} is used to represent
     mortality caused by background sources such as illness or age. The
     rate of senescence mortality (in 1/year) is given by:

     \deqn{\mu_{sen.i}(w) = k_{sen}\left(
                                 \frac{log_{10}(w)}{log_{10}(w_{max.i})}
                                 \right)^{p_{sen}}}
          {\mu_{sen.i}(w) = k_{sen}
          (log_{10}(w)/log_{10}(w_{max.i}))^{p_{sen}}}

     where \eqn{k_{sen}} is the rate of senescence mortality, \eqn{p_sen}
     defines the slope of the senescence curve, \eqn{w_max.i} maximum body
     size of group \eqn{i} in grams.
