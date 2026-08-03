# Mass-specific algae consumption rate

This mass-specific consumption rate is used in
[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md)
to calculate the algae biomass at the next time step. To get the
non-mass-specific consumption rate, use
[`getAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeConsumption.md).

## Usage

``` r
algae_consumption(params, n = params@initial_n, rates = getRates(params))
```

## Arguments

- params:

  MizerParams

- n:

  A matrix of current species abundances (species x size)

- rates:

  A list of rates as returned by
  [`getRates()`](https://sizespectrum.org/mizer/reference/getRates.html)

## Value

The mass-specific consumption rate of algae in grams per year.

## Details

The rho parameter for herbivorous fish groups is stored in
`other_params(params)$algae$rho`

## Algae consumption

This rate deliberately does not depend on feeding level or the
`satiation` species parameter (contrast with
[`detritus_consumption()`](https://cmbeese.github.io/mizerReef/reference/detritus_consumption.md),
which does): in mizerReef, satiation-mediated consumption is exclusive
to detritivory. Increases in herbivorous fish density following coral
bleaching events suggest that reef herbivores respond to increased food
availability without regulating their consumption (Ledlie et al. 2007;
Pratchett et al. 2008; Khalil et al. 2013; Elma et al. 2023), and
Caribbean herbivores have been observed to fill their gut up to three
times a day (Ferreira et al. 1998; Kopp et al. 2010). Algal depletion is
therefore modelled as driven by continuous grazing pressure rather than
by any individual consumer's satiation state.
[`getAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeConsumption.md)
reports the feeding-level-adjusted rate actually ingested by each
species for diagnostic purposes, but that adjusted rate is not what
depletes the algae pool or what
[`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
use for tuning.

The rate at which herbivorous consumer groups encounter algae biomass
\\E\_{i.A}(w)\\ is controlled by the parameter \\\rho\_{A.i}\\. It
scales with the size of the consumer raised to an allometric exponent
\\m\_{alg}\\ which is taken from empirical data.

\$\$E\_{i.A}(w)=\rho\_{i.A}\\ w^{m\_{alg}}\\B_A\$\$

The mass specific consumption rate then accounts for the preference of
functional group \$i\$ for algae, \\\theta\_{i.A}\\. This gives the
mass-specific algae consumption rate:

\$\$c_A = \sum_i\int\rho\_{i.A}\\ w^{m\_{alg}}
N_i(w)\theta\_{i.A}\\dw\$\$

## Examples

``` r
data(caribbean_3_model)
algae_consumption(caribbean_3_model)
#> [1] 9.21051e+12
```
