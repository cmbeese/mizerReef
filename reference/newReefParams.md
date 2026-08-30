# Set up parameters for a mizerReef model

Sets up a multi-species size spectrum model with additional unstructured
resource components, senescence mortality, and predation refuge.

## Usage

``` r
newReefParams(
  species_params,
  interaction = NULL,
  crit_feed = NULL,
  min_w_pp = NA,
  w_pp_cutoff = 1,
  n = 0.75,
  new_refuge = FALSE,
  method,
  method_params,
  refuge_user = NULL,
  blocked_pred = NULL,
  satiation = NULL,
  a_bar = NULL,
  b_bar = NULL,
  w_settle = NULL,
  max_protect = NULL,
  tau = NULL,
  use_dummy_fish_bins = TRUE,
  degrade = FALSE,
  bleach_time = 2,
  trajectory = NULL,
  deg_scale = 1,
  algae_boost = FALSE,
  algae_growth_boost = c(1.11, 1.11, 1.11, 1.11),
  algae_capacity_boost = c(2),
  UR_interaction = NULL,
  use_UR_cc = FALSE,
  initial_algae_growth = NULL,
  algae_capacity = NULL,
  detritus_capacity = NULL,
  sen_decomp = NULL,
  ext_decomp = NULL,
  external = NULL,
  algae_colour = "darkseagreen3",
  detritus_colour = "plum4",
  resource_color = "lightseagreen",
  ext_mort_params = NULL,
  include_ext_mort = TRUE,
  include_sen_mort = TRUE,
  z0pre = 0.2,
  info_level = mizer::default_info_level(),
  ...
)
```

## Arguments

- species_params:

  A species parameter data frame containing at least the name of each
  species, their observed abundances, and their maximum size.

- interaction:

  The species specific interaction matrix, \\\theta\_{ij}\\

- crit_feed:

  Critical feeding level

- min_w_pp:

  Minimum size of plankton in grams

- w_pp_cutoff:

  Maximum size of plankton in grams, default to 1 g

- n:

  Allometric growth exponent (also used as metabolic exponent p)

- new_refuge:

  Logical. If TRUE, indicates this refuge profile is being set for use
  in a new simulation (prevents algae/detritus from being re-tuned
  during
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)).
  Default FALSE.

- method:

  Character. The method for setting up benthic refuge. One of
  "sigmoidal", "binned", "competitive", or "noncomplex". See Details.
  **Required.**

- method_params:

  Data frame or named list. Specifies parameters required for the chosen
  method:

  - For "sigmoidal": must include `L_refuge` (numeric, length at which
    refuge becomes scarce, cm; **no default**) and `prop_protect`
    (numeric, max proportion protected, default: 0.98).

  - For "binned": must include `start_L` (numeric, start length, cm;
    **no default**), `end_L` (numeric, end length, cm; **no default**),
    and `prop_protect` (numeric, proportion protected, default: 0.98).

  - For "competitive": must include `start_L`, `end_L` (as above), and
    `refuge_density` (numeric, refuges per size bin, no/m^2; **no
    default**).

  - For "noncomplex": no parameters required.

- refuge_user:

  Logical vector (length = number of species). Indicates which groups
  use refuge. If not present in `species_params`, must be provided.
  Defaults to FALSE.

- blocked_pred:

  Optional. Logical vector (length = number of species). Indicates
  whether the predator is blocked by refuge for this species. TRUE means
  hunting is blocked by refuge; FALSE means the species can access prey
  within refuge (e.g. eels). Defaults to FALSE.

- satiation:

  Logical vector (length = number of species). Indicates which groups
  are subject to satiation. In mizerReef, satiation is intended to be
  exclusive to detritivory (see
  [`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)'s
  "Algae consumption" section) – if not provided, defaults are set
  automatically: TRUE only for detritivores (species with positive
  `interaction_detritus` that do NOT also graze algae, i.e.
  `interaction_algae` is zero or absent, and that do not eat other
  species, i.e. row sum of interaction matrix is 0); FALSE for every
  other species, including carnivores, pure algae/plankton grazers, and
  species that consume both algae and detritus (their diet is
  herbivore-like, so they default to the unregulated, herbivore-style
  behaviour rather than the detritivore-style one). A warning is issued
  if defaults are used.

- a_bar:

  Numeric. Length-weight conversion parameter for dummy fish.
  **Default:** 0.025. If any species is missing an 'a' parameter, the
  value of a_bar is used for that species and a warning is issued.

- b_bar:

  Numeric. Length-weight exponent for dummy fish. **Default:** 3. If any
  species is missing a 'b' parameter, the value of b_bar is used for
  that species and a warning is issued.

- w_settle:

  Numeric. Minimum weight of fish protected by refuges at measured scale
  (grams). **Default:** 0.1.

- max_protect:

  Numeric. Maximum proportion of fish protected by refuge (0–1).
  **Default:** 0.98.

- tau:

  Numeric. Proportion of fish with access to refuge expected to utilize
  it (0–1). **Default:** 1.

- use_dummy_fish_bins:

  Logical. Controls how refuge bin boundaries and thresholds are
  calculated for sigmoidal, binned, and competitive methods:

  - TRUE (default, legacy behavior): Use dummy fish parameters (`a_bar`,
    `b_bar`) to set bin boundaries/thresholds.

  - FALSE: Use each species' own length-weight parameters (`a`, `b`).
    Set according to your data collection method. The setting is stored
    in `params@other_params$refuge_params$use_dummy_fish_bins` and used
    by
    [`getRefuge()`](https://cmbeese.github.io/mizerReef/reference/getRefuge.md).

- degrade:

  Logical. Whether to enable habitat degradation during projections.
  Default FALSE. See
  [`setDegradation()`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md)
  for details on degradation parameters.

- bleach_time:

  The year of the simulation to implement bleaching. Defaults to year 2.

- trajectory:

  Optional character string naming the degradation scenario. Used for
  documentation and identification. No functional effect on the
  simulation.

- deg_scale:

  A matrix (refuge bins x time) giving scaling factors for refuge
  density. Column 1 represents the bleaching year (initial impact), and
  subsequent columns represent years 1, 2, 3... post-bleaching. Values
  are multiplied by the previous timestep's refuge density. Default
  scaling matrices for 15 years with "rubble", "algae", and "recovery"
  trajectories are included as data objects in the package.

- algae_boost:

  Logical. Should algae growth and/or carrying capacity be adjusted in
  response to bleaching? Default is FALSE.

- algae_growth_boost:

  Numeric vector. Multipliers for algae growth rate at bleaching
  (element 1) and during post-bleaching years (subsequent elements).
  Only used if algae_boost = TRUE. Vector length determines duration:
  length 4 means bleaching year plus 3 post-bleaching years. Default is
  c(1.11, 1.11, 1.11, 1.11) for 4 years of boosting.

- algae_capacity_boost:

  Numeric vector. Multipliers for algae carrying capacity at bleaching
  (element 1) and during post-bleaching years (subsequent elements).
  Only used if algae_boost = TRUE. If shorter than algae_growth_boost,
  will be padded with 1s (no change). Default is c(2.0) (only boosts at
  bleaching year).

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

- initial_algae_growth:

  Numeric. The fixed, literature-informed growth rate of algae in
  grams/m^2/year, passed on to
  [`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)
  as `algae_growth_initial` (see there for the literature basis for the
  default of 2e3). Held constant rather than retuned to match
  consumption – see
  [`?setAlgaeParams`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)'s
  "Algae as an unstructured resource" section.

- algae_capacity:

  Numeric. Carrying capacity for algae biomass in grams per year.
  Default is 1.

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

- algae_colour:

  Character. Colour to use for algae in plots. Default is
  "darkseagreen3".

- detritus_colour:

  Character. Colour to use for detritus in plots. Default is "plum4".

- resource_color:

  Character. Colour to use for the resource line in plots. Default is
  "lightseagreen".

- ext_mort_params:

  Optional. A named list or matrix with columns/names 'nat_mort',
  'sen_prop', and 'sen_curve'. If NULL, defaults are used. All values
  must be numeric and non-negative.

- include_ext_mort:

  A boolean value that indicates whether the user wants to use default
  external mortality. Defaults to TRUE.

- include_sen_mort:

  A boolean value that indicates whether the user wants to use default
  senescence mortality. Defaults to TRUE.

- z0pre:

  If `include_ext_mort`is FALSE, the external mortality rate for each
  species calculated as z0pre \* w_max ^ z0exp. z0exp defaults to 1-n
  where n is the given allometric scaling exponent and z0pre defaults to
  0.2.

- info_level:

  How much mizer should say about the choices it makes here. Level 1
  keeps only the reports that tell you something went differently from
  how you asked; 0 is silence. See
  [`mizer::default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.html).

- ...:

  Extra parameters to be passed to
  [`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.html)

## Value

An object of type
[MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)

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

## Setting the refuge profile

     The mizerReef package provides three methods to define the refuge profile.

     \itemize{

     \item **Sigmoidal Method**: \cr

         This method is preferred for data-poor reefs or reefs where the refuge
         distribution is unknown. It is also ideal for systems where only one
         species is expected to be utilizing refuge. The sigmoidal method defines
         a smooth transition in refuge availability around a threshold body size.

         The threshold for refuge can be set in two ways, depending on how your
         refuge data was collected:

         - If `use_dummy_fish_bins = FALSE`, the threshold weight is calculated as:
                 \deqn{ W_{i.refuge} = a_i \cdot L_{refuge}^{b_i} }
          where \eqn{a_i} and \eqn{b_i} are the length-weight parameters for species i.

         - If `use_dummy_fish_bins = TRUE`, the weight threshold is:

                 \deqn{ W_{refuge} = a_{bar} \cdot L_{refuge}^{b_{bar}} }

          The proportion of fish with access to refuge is then given by:
                 \deqn{ R_{j}(w_p) = \frac{r}{1 + e^{\Delta(w - W_{refuge})}} }
         where $r$ is the maximum proportion protected, $\Delta$ is the slope,
         $w$ is body weight, and $W_{refuge}$ is the threshold (species-specific or dummy fish)

         For this method, `method_params` should contain columns named
         `prop_protect` and `L_refuge` that give the values for \eqn{r}
         and the length at which refuge becomes scarce in cm.

     \item **Binned Method**: \cr

         This method is appropriate for theoretical applications
         and does not rely on empirical data. It sets refuge to a constant
         proportion of fish within a given size range. The proportion of fish
         in group \eqn{j} with access to refuge is given by

         \deqn{ R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~] }{
                  R_j(w_p) = r_k ~~~~~~~ w_p ∈ (~w_{k-1}, w_k~] }

         where \eqn{r_k} is the proportion of fish with access to refuge in
         size class \eqn{k}.

         For this method, `method_params` should contain columns named
         `start_L` and `end_L` which contain the starting and ending lengths [cm]
         of each size bin and `prop_protect`, the proportion of fish protected
         within each corresponding size bin.

     \item **Competitive Method**: \cr
         This method is appropriate when refuge density data is available for
         the modelled reef. The refuge density describes the distribution of
         refuges \eqn{(no./m^2)} across predefined fish body size categories.
         The proportion of fish in size class \eqn{k} with access to refuge
         is given by

         \deqn{R_{j}(w_p) = \tau \cdot \frac{ \eta_{k} }
                                                { \sum_i \int_{w_{k-1}}^{w-k} N_i(w) \, dw}}{
                  R_{j}(w_p) = \tau \eta_{k} /
                              ( \sum_i \int_{w_{k-1}}^{w-k} N_i(w) \, dw ) }

         where \eqn{ \tau } is the proportion of fish with access to refuge that
         are expected to actually utilize  it, \eqn{ \eta_{k}} is the density of
         refuges in size range \eqn{(w_{k-1}, w_k]} and
         \eqn{\sum_{i} \int_{w_{k-1}}^{w_k} N_i(w)~dw} gives the density
         of fish from any group in size range \eqn{(w_{k-1}, w_k]}.
         This represents the density of competitors for refuges in
         size class \eqn{k}.

         For this method, `method_params` should contain columns named
         `start_L`and `end_L` which contain the starting and ending lengths [cm]
         of each size bin and `refuge_density`, the number of refuges available
         in each size bin (no/m^2).

     }

Users can also set a noncomplex reef with no habitat refuge. This option
is convenient for finding steady state parameters and is the default
when no parameters are provided.

This function checks that the supplied refuge parameters are valid, adds
relevant columns to the `species_params` data frame, and stores refuge
parameters in the `refuge_params` slot of the `params` object.

Refuge profile parameters can be input in a spreadsheet program and
saved as a .csv file. The data can then be read into R using the command
[`read.csv()`](https://rdrr.io/r/utils/read.table.html).

## Examples

``` r
params <- newReefParams(
    species_params = caribbean_3_species,
    interaction = caribbean_3_interaction,
    method = "binned",
    method_params = tuning_profile
)
#> ℹ No h provided for some species, so using age at maturity to calculate it.
#> ℹ For species where no growth information is available the parameter h has been set to h = 30.
#> ℹ Using z0 = z0pre * w_inf ^ z0exp for calculated z0 values.
#> ℹ Using f0, h, lambda, kappa and the predation kernel to calculate gamma.
class(params)
#> [1] "mizerReef"
#> attr(,"package")
#> [1] ".GlobalEnv"
```
