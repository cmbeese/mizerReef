# Set the refuge profile parameters

Sets and validates refuge parameters in the mizerParams object, allowing
flexible calibration to match ecological data and measurement
approaches. Supports multiple methods for defining refuge profiles and
options for calibrating bin boundaries using either species-specific or
dummy fish length-weight parameters.

## Usage

``` r
setRefuge(
  params,
  method,
  method_params = NULL,
  refuge_user = NULL,
  blocked_pred = NULL,
  satiation = NULL,
  a_bar = NULL,
  b_bar = NULL,
  w_settle = NULL,
  max_protect = NULL,
  tau = NULL,
  use_dummy_fish_bins = TRUE,
  ...
)
```

## Arguments

- params:

  MizerParams object

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

- ...:

  Unused.

## Value

A MizerParams object with updated refuge parameters

## Details

Refuge profiles account for the protective behavior of prey living in
high-complexity environments (e.g. coral reefs) with access to predation
refuge. The refuge profile defines the proportion of fish within
user-defined length bins that are protected from being encountered by a
predator.

A unique refuge profile is generated for each predator group × prey
group × prey size combination based on the given refuge profile
parameters and four values from `params@species_params`: length-weight
conversion values `a` and `b`, `refuge_user` (TRUE for groups that
utilize predation refuge), and `blocked_pred` (FALSE for predator groups
whose body shape or predatory strategy allow them to access fish within
refuge, e.g. eels).

The maximum proportion of fish protected by refuge in any size class is
set by `max_protect` to ensure some food is always available to
predators.

The refuge profile is used when calculating the food encounter rate in
[`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md)
and the predation mortality rate in
[`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md).
Its entries are dimensionless values between 0 and 1, representing the
proportion of fish in the corresponding prey and size categories that
are hidden within refuge and thus cannot be encountered by predators. If
no refuge is available, predator-prey interactions are determined
entirely by size-preference.

## Defining refuge threshholds

     Defining the refuge profile for a given system is nontrivial. There are many
     ways to measure and define predation refuge on reefs. The way you parameterise
     the model should reflect your data collection method and ecological context.

     The `use_dummy_fish_bins` argument determines how refuge bin boundaries are
     set:

     - If TRUE (default), bin boundaries are determined by weight, calculated
       using dummy fish length-weight parameters (`a_bar`, `b_bar`). All species
       share the same weight boundaries for bins, but the lengths of fish in each
       bin are calculated for each species using their own length-weight
       parameters (`a`, `b`). This means fish in the same bin may be of very
       different lengths depending on species. This is appropriate when refuge
       capacity is measured using dummy fish of known size and is the
       preferred method.

     - If FALSE, bin boundaries are determined by length, so all species in a bin
       have the same length, but their weights are calculated using their own
       length-weight parameters (`a`, `b`). This means fish in the same bin may
       have very different weights depending on species. This is appropriate when
       refuge hole entrances are measured for each species.

This choice applies to all refuge methods (sigmoidal, binned,
competitive) and is stored in
`params@other_params$refuge_params$use_dummy_fish_bins`. Select the
option that matches your measurement approach and ecological realism.

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

## See also

[`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md),
[`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md),
[`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md),
[`setDetritusParams()`](https://cmbeese.github.io/mizerReef/reference/setDetritusParams.md)

## Examples

``` r
data(caribbean_3_model)
data(tuning_profile)
params <- setRefuge(caribbean_3_model,
    method = "binned", method_params = tuning_profile
)
```
