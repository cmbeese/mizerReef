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

  A list of rates as returned by `getRates()`

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

\$\$p\_{D.f} = \sum_i(1-\alpha_i)\int E_i(w)\\dw\$\$

\\p\_{D.d}\\ comes from the biomass of fish that die as a result of
external mortality. External mortality includes local deaths that lead
to detritus but also deaths due to predation by species that are not
explicitly modelled, for example transient predators, mammals, or sea
birds. Thus, only a proportion `prop_decomp` of this material decomposes
to detritus. The detritus production from decomposing dead organisms is
given by:

\$\$p\_{D.d} = \sum_i\int\mu\_{seni.i}(w)N_i(w)w\\dw +
\mathtt{prop\\decomp}\\ \sum_i\int\mu\_{nat.i}(w)N_i(w)w\\dw\$\$

\\p\_{D.ext}\\ is the rate at which detritus enters the system from
unmodelled or external sources. For coral reefs, this includes detritus
produced by sponges and coral mucous as well as waste material that
sinks in from the pelagic zone. This rate is a model parameter
independent of any other model component. It is set so that production
and consumption are equal for the chosen steady state abundances.
