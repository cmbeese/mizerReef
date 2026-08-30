# Change the refuge parameters for a model

This is a wrapper function for the
[`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)
and
[`getRefuge()`](https://cmbeese.github.io/mizerReef/reference/getRefuge.md)
functions that allows users to easily change refuge parameters on an
existing mizer model. This will take it out of steady state, and users
should run
[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
after this function to return to steady state.

## Usage

``` r
newRefuge(
  params,
  new_refuge = FALSE,
  new_method = NULL,
  new_method_params = NULL,
  new_L_refuge = NULL,
  new_prop_protect = NULL,
  scale_bin = NULL,
  info_level = mizer::default_info_level(),
  ...
)
```

## Arguments

- params:

  a mizer object

- new_refuge:

  A boolean value that states whether this new refuge profile is being
  used for simulation. Determines whether algae and detritus production
  are tuned when the model is run to steady state. Defaults to FALSE.
  Set this to `TRUE` when comparing fish outcomes under two refuge
  profiles (e.g. before/after a habitat complexity or MPA-style change)
  and you want to isolate the effect of the refuge change on fish by
  holding the algae/detritus resource base fixed at its current values,
  rather than letting it silently readjust and mask the comparison.
  Leave at the default `FALSE` for ordinary steady-state tuning, where
  algae/detritus should retune to match the (possibly still-changing)
  fish community. This rationale is inferred from the code's behaviour,
  not confirmed against original design notes – verify it matches your
  intended use before relying on it for a scenario comparison.

- new_method:

  The new method to be used for setting the refuge profile. Options are
  "sigmoidal", "binned", "competitive", or "noncomplex". If no method is
  provided, this defaults to the same method as is currently being used
  in the simulation.

- new_method_params:

  A data frame or list specifying parameters required for the chosen
  method. For 'sigmoidal', must include 'L_refuge' and 'prop_protect'.
  For 'binned' and 'competitive', must include 'start_L', 'end_L', and
  'prop_protect' (or 'refuge_density' for competitive). Only necessary
  if changing methods.

- new_L_refuge:

  To be used with "sigmoidal" method only. The new value for the length
  at which refuge becomes scarce in cm.

- new_prop_protect:

  To be used with "sigmoidal" method only. The new value for the maximum
  proportion of fish to protect.

- scale_bin:

  To be used with the "binned" or "competitive" method only. A number or
  vector (length = number of bins) to multiply the `prop_protect`
  (binned) or `refuge_density` (competitive) values for each size bin.
  Changes the proportion of fish protected. If a single value is given,
  it is applied to all bins.

- info_level:

  How much mizer should say about the choices it makes here. Level 1
  keeps only the reports that tell you something went differently from
  how you asked; 0 is silence. See
  [`mizer::default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.html).

- ...:

  Unused

## Value

A mizer object with updated refuge profiles

## Details

At least one of the arguments new_method, new_method_params,
new_L_refuge, new_prop_protect, or scale_bin must be provided. For
'binned' and 'competitive' methods, scale_bin must be a value or vector
matching the number of bins.

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
