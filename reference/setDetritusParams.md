# Set detritus parameters for mizerReef

Set detritus parameters for mizerReef

## Usage

``` r
setDetritusParams(
  params,
  detritus_capacity = NULL,
  sen_decomp = NULL,
  ext_decomp = NULL,
  external = NULL,
  UR_interaction = NULL,
  use_UR_cc = FALSE,
  detritus_colour = "plum4"
)
```

## Arguments

- params:

  A `MizerParams` object.

- detritus_capacity:

  Numeric. Carrying capacity for detritus biomass in grams per year.
  Default is 1.

- sen_decomp:

  Numeric. Proportion of decomposing mass from senescence mortality that
  becomes detritus. Default is 0.8.

- ext_decomp:

  Numeric. Proportion of decomposing mass from external mortality that
  becomes detritus. Default is 0.2.

- external:

  Numeric. Rate at which detritus biomass sinks from the pelagic zone
  (grams per year). Default is 1.

- UR_interaction:

  Optional. A named list or array with one or more resource interaction
  vectors (e.g. interaction_algae, interaction_detritus,
  interaction_sponge), each of length equal to the number of species. If
  NULL, will use columns in species_params or set to zero. All values
  must be numeric and between 0 and 1.

- use_UR_cc:

  Logical. Whether to implement a carrying capacity for detritus.
  Default is FALSE. This flag is stored in the other_params slot.

- detritus_colour:

  Character. Colour to use for detritus in plots. Default is "plum4".

## Value

A `MizerParams` object with updated detritus parameters (in
`other_params$detritus`).

## Details

All detritus-related parameters (capacity, decomposition rates, external
input) are stored in the `detritus` component of `other_params` (i.e.
`other_params(params)$detritus`), the same location
[`mizer::getComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)/[`mizer::removeComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)
use. Resource interaction strengths are set in the species_params data
frame. This function supports flexible multi-resource interaction via
the UR_interaction argument.

## Detritus as an unstructured resource

     Detritus in mizerReef is modeled as a non-size-structured resource
     produced by the decomposition of organic materials. Detritus is consumed
     by detritivores and benthic invertebrates. This function sets the
     carrying capacity, decomposition proportions, and interaction strengths
     for detritus, supporting flexible diet preferences.

     The interaction strength (\eqn{\theta_{i,detritus}}) for each species \eqn{i}
     determines how strongly that group feeds on detritus. This can be set via the
     `interaction_detritus` column in the species parameter data frame, or directly
     via the `UR_interaction` argument. If neither is provided, all interaction
     strengths are set to zero and a warning is issued.

     The carrying capacity (`detritus_capacity`) limits the maximum standing
     stock of detritus.

     The proportions of decomposing mass from senescence (`sen_decomp`) and
     external mortality (`ext_decomp`) that become detritus are set here,
     with typical defaults of 0.8 and 0.2, respectively. The rate at which
     detritus sinks from the pelagic zone (`d_external_initial`) is also set,
     controlling detritus input from external sources. These values may be
     reset by [reefSteady()] to ensure steady state abundances match observed
     or target values.

     Carrying capacity can be toggled with `use_UR_cc`. When enabled, detritus
     biomass will be limited by the specified capacity.

     Note: Interaction with size-structured resources, such as plankton, is
     set with the resource_interaction column of the species parameters dataframe.

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
(external mortality, default 0.2; see `setDetritusParams()`). The
detritus production from decomposing dead organisms is given by:

\$\$p\_{D.d} = \mathtt{sen\\decomp}\\
\sum_i\int\mu\_{seni.i}(w)N_i(w)w\\dw + \mathtt{ext\\decomp}\\
\sum_i\int\mu\_{nat.i}(w)N_i(w)w\\dw\$\$

\\p\_{D.ext}\\ is the rate at which detritus enters the system from
unmodelled or external sources. For coral reefs, this includes detritus
produced by sponges and coral mucous as well as waste material that
sinks in from the pelagic zone. This rate is a model parameter
independent of any other model component. It is set so that production
and consumption are equal for the chosen steady state abundances.

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

## Examples

``` r
data(caribbean_3_model)
params <- setDetritusParams(caribbean_3_model, external = 10)
```
