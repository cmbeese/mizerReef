# Set algae parameters for mizerReef

Set algae parameters for mizerReef

## Usage

``` r
setAlgaeParams(
  params,
  algae_growth_initial = NULL,
  algae_capacity = NULL,
  UR_interaction = NULL,
  use_UR_cc = FALSE,
  algae_colour = "darkseagreen3",
  info_level = mizer::default_info_level()
)
```

## Arguments

- params:

  A `MizerParams` object.

- algae_growth_initial:

  Numeric. The fixed, literature-informed growth rate of algae in
  grams/m^2/year, held constant (not retuned to match consumption – see
  the "Algae as an unstructured resource" section above). Default is
  2e3, near the upper end of the range implied by converting Caribbean
  algal turf net production estimates of ~110-475 g dry weight/m^2/year
  (Carpenter, R.C. (1986). Partitioning herbivory and its effects on
  coral reef algal communities. Ecological Monographs, 56, 345-363;
  consistent with the meta-analysis in Tebbett, S.B. & Bellwood, D.R.
  (2021). Algal turf productivity on coral reefs: A meta-analysis.
  Marine Environmental Research, 168, 105311) to the wet-mass units used
  elsewhere in mizerReef, using an approximate 20% dry-matter content
  for turf algae.

- algae_capacity:

  Numeric. Carrying capacity for algae biomass in grams per year.
  Default is 1.

- UR_interaction:

  Optional. A named list or array with one or more resource interaction
  vectors (e.g. interaction_algae, interaction_detritus,
  interaction_sponge), each of length equal to the number of species. If
  NULL, will use columns in species_params or set to zero. All values
  must be numeric and between 0 and 1.

- use_UR_cc:

  Logical. Whether to implement a carrying capacity for all unstructured
  resources. Default is FALSE. This flag is stored in the other_params
  slot.

- algae_colour:

  Character. Colour to use for algae in plots. Default is
  "darkseagreen3".

- info_level:

  How much mizer should say about the choices it makes here. Level 1
  keeps only the reports that tell you something went differently from
  how you asked; 0 is silence. See
  [`mizer::default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.html).

## Value

A `MizerParams` object with updated algae parameters (in
`other_params$algae`).

## Details

All algae-related parameters (growth rate, capacity) are stored in the
`algae` component of `other_params` (i.e. `other_params(params)$algae`),
the same location
[`mizer::getComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)/[`mizer::removeComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)
use. Resource interaction strengths are set in the species_params data
frame. This function supports flexible multi-resource interaction via
the UR_interaction argument.

## Algae as an unstructured resource

     mizerReef supports algae as a non-size-structured resource,
     consumed primarily by herbivorous fish. This function sets
     the initial growth rate, system carrying capacity, and
     interaction strengths for algae, allowing for flexible
     diet preferences.

     The interaction strength (\eqn{\theta_{i,algae}}) for each
     species \eqn{i} determines how strongly that group feeds on
     algae. This can be set via the `interaction_algae` column in
     the species parameter data frame, or directly via the
     `UR_interaction` argument. If neither is provided, all
     interaction strengths are set to zero and a warning is issued.

     The initial growth rate (`algae_growth_initial`) and carrying
     capacity (`algae_capacity`) control the baseline production
     and maximum standing stock of algae, respectively. Unlike
     detritus, algal production on a reef is real primary production
     and is not driven by grazer demand, so `algae_growth_initial` is
     treated as a fixed, literature-informed constant: it is *not*
     reset by [reefSteady()]. Instead [reefSteady()] (via [tuneUR()]/
     [tuneUR_cc()]) solves for the algae *biomass* that is at steady
     state for this fixed production rate and the model's current
     abundances, so that, all else equal, reducing grazing pressure
     increases the resulting algae biomass rather than reducing
     production to compensate.

     Carrying capacity can be toggled with `use_UR_cc`. When enabled,
     algae biomass will be limited by the specified capacity.

     Note: Interaction with size-structured resources, such as plankton,
     is set with the resource_interaction column of the species parameters
     dataframe.

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
params <- setAlgaeParams(caribbean_3_model, algae_growth_initial = 2500)
```
