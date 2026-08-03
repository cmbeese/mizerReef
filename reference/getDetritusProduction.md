# Detritus production rate

Gives a named vector with the rates at which different components of the
ecosystem produce detritus:

## Usage

``` r
getDetritusProduction(params, n = params@initial_n, rates = getRates(params))
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

A vector with named entries "feces", "decomp", and "external", giving
the rates at which detritus biomass is produced by each of these sources
in grams per year.

## Details

1.  consumed biomass not assimilated by predators ("feces"),

2.  decomposing dead organisms ("decomp"),

3.  the pelagic zone ("external").

This function returns a vector with the individual contributions for
each source. These can be summed with
[`sum()`](https://rdrr.io/r/base/sum.html) to get the total detritus
production rate.

## Detritus production

The rate \\p_D\\ at which detritus biomass is produced by the ecosystem
has contributions from three sources:

\$\$p_D = p\_{D.f} + p\_{D.d} + p\_{D.ext}\$\$

\\p\_{D.f}\\ comes from the biomass that is consumed but not assimilated
and is given by:

\$\$p\_{D.f} = \sum_i(1-\alpha_i)\int (1-f_i(w))\\E_i(w)\\N_i(w)\\dw\$\$

where \\f_i(w)\\ is the feeding level (see
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)'s
"Algae consumption" section for how `satiation` controls it), so that
\\(1-f_i(w))\\E_i(w)\\ is the biomass actually consumed (as opposed to
merely encountered) – unlike algae consumption, which deliberately
ignores feeding level (see
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)),
detritus's egestion term uses the same feeding-level-adjusted
consumption rate as
[`getDetritusConsumption()`](https://cmbeese.github.io/mizerReef/reference/getDetritusConsumption.md)
and
[`detritus_consumption()`](https://cmbeese.github.io/mizerReef/reference/detritus_consumption.md).

\\p\_{D.d}\\ comes from the biomass of fish that die, combining two
mortality sources that each decompose to detritus at their own rate:
senescence mortality (see
[`getSenMort()`](https://cmbeese.github.io/mizerReef/reference/getSenMort.md))
and external mortality, i.e. local deaths that lead directly to detritus
as well as deaths due to predation by species that are not explicitly
modelled, for example transient predators, mammals, or sea birds. Only a
proportion of each source's dead biomass decomposes to detritus, set
independently by `sen_decomp` (senescence, default 0.8) and `ext_decomp`
(external mortality, default 0.2; see
[`setDetritusParams()`](https://cmbeese.github.io/mizerReef/reference/setDetritusParams.md)).
The detritus production from decomposing dead organisms is given by:

\$\$p\_{D.d} = \mathtt{sen\\decomp}\\
\sum_i\int\mu\_{seni.i}(w)N_i(w)w\\dw + \mathtt{ext\\decomp}\\
\sum_i\int\mu\_{nat.i}(w)N_i(w)w\\dw\$\$

\\p\_{D.ext}\\ is the rate at which detritus enters the system from
unmodelled or external sources. For coral reefs, this includes detritus
produced by sponges and coral mucous as well as waste material that
sinks in from the pelagic zone. This rate is a model parameter
independent of any other model component. It is set so that production
and consumption are equal for the chosen steady state abundances.

## Examples

``` r
data(caribbean_3_model)
getDetritusProduction(caribbean_3_model)
#>     feces    decomp  external 
#>  176.0479  119.4528 -110.7587 
```
