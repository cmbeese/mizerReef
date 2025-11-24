# Checks unstructured resource parameters and interaction matrix

Checks unstructured resource parameters and interaction matrix

## Usage

``` r
setURParams(
  params,
  UR_interaction = NULL,
  initial_algae_growth = NULL,
  carry_capacity = FALSE,
  algae_capacity = NULL,
  detritus_capacity = NULL,
  sen_decomp = NULL,
  ext_decomp = NULL,
  initial_d_external = NULL
)
```

## Arguments

- params:

  MizerParams object

- UR_interaction:

  Interaction matrix for unstructured resources (species x resource)

- initial_algae_growth:

  The initial growth rate of algae in grams/m^2/year. This value is
  reset to match consumption in the
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  function so that steady state abundances match given values.

- carry_capacity:

  A boolean value that indicates whether the user wants to implement a
  carrying capacity for unstructured resources. Default is FALSE

- algae_capacity:

  The carrying capacity of the system for algae biomass in grams per
  year.

- detritus_capacity:

  The carrying capacity of the system for detritus biomass in grams per
  year.

- sen_decomp:

  The proportion of decomposing mass from senescence mortality that
  decomposes to become part of the detritus pool. Defaults to 0.8.

- ext_decomp:

  The proportion of decomposing mass from external mortality that
  decomposes to become part of the detritus pool. Defaults to 0.2.

- initial_d_external:

  The rate at which detritus biomass sinks from the pelagic zone and
  becomes part of the detritus pool in grams per year. This value is
  reset to make up any differences in consumption and production in the
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  function so that steady state abundances match observed values.

## Value

`setURParams` MizerParams object with updated unstructured resource
parameters

## Adding unstructured resources

     mizerReef supports two resource spectra that are not size- structured.
     Algae are consumed by herbivorous fish, while detritus is consumed by
     herbivorous fish and benthic invertebrates. This function sets the
     interaction matrix for these resources as well as any default
     parameters necessary to structure them.

     The resource interaction matrix \eqn{\theta_{ki}} modifies the
     interaction of each functional group \eqn{i} with each unstructured
     resource \eqn{k} in the model. This can be used for example to allow
     for different diet preferences on each unstructured resource.

     Note that interaction with size structured resources, such as
     plankton, is still set with the resource_interaction column of
     the species parameters dataframe.

## Detritus consumption

The rate at which detritivorous consumer groups encounter detrital
biomass \\E\_{i.D}(w)\\ is controlled by the parameter \\\rho\_{D.i}\\.
It scales with the size of the consumer raised to an allometric exponent
\\m\_{det}\\ which is taken to be the same as the scaling exponent of
the maximum intake rate for fish consumers.

\$\$E\_{i.D}(w)=\rho\_{i.D}\\ w^{m\_{det}}\\B_D \$\$

The mass specific consumption rate then accounts for the preference of
functional group \$i\$ for detritus, \\\theta\_{i.D}\\ and the feeding
level \\f_i(w)\\. This gives the mass-specific detritus consumption
rate:

\$\$c_D = \sum_i\int\rho\_{i.D}\\ w^{m\_{det}} N_i(w)
\left(1-f_i(w)\right) \theta\_{i.D}\\dw\$\$

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

## Algae consumption

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
