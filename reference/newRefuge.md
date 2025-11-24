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
  ...
)
```

## Arguments

- params:

  a mizer object

- new_refuge:

  A boolean value that states whether this new refuge profile is being
  used for simulation. Determines whether algae and detritus production
  and tuned when the model is run to steady state. Defaults to false.

- new_method:

  The new method to be used for setting the refuge profile. Options are
  "sigmoidal", "binned", "competitive", or "noncomplex". If no method is
  provided, this defaults to the same method as is currently being used
  in the simulation.

- new_method_params:

  A data frame containing values specific to each method for calculating
  refuge. Only necessary if changing methods.

- new_L_refuge:

  To be used with "sigmoidal" method only. The new value for the length
  at which refuge becomes scarce in cm.

- new_prop_protect:

  To be used with "sigmoidal" method only. The new value for the maximum
  proportion of fish to protect.

- scale_bin:

  To be use with the "binned" method only. A number or vector of numbers
  to multiply the `prop_protect` values by for each size bin. Changes
  the proportion of fish protected.

- ...:

  Unused

## Value

A mizer object with updated refuge profiles

## Setting the refuge profile

Refuge profiles account for the protective behavior of prey living in
high complexity environments (e.g. coral reefs) with access to predation
refuge. The refuge profile defines the proportion of fish within
user-defined length bins that are protected from being encountered by a
predator.

mizerReef determines how much fish weight a refuge can hold by
converting user provided fish length bins to weight bins with average
length to weight conversion parameters \\\bar{a}\\ and \\\bar{b}\\,
which describe the length to weight relationship for the fish dummies
used during data collection. The default values are \\\bar{a} = 0.025,
\bar{b} = 3\\. Assuming that fish in shelters are neutrally buoyant, the
mass of the largest fish in each size category is proportional to the
volume of refuges classified in that category.

A unique refuge profile is generated for each predator group x prey
group x prey size combination based on the given refuge profile
parameters as well as four values from `params@species_params`: length
to weight conversion values `a` and `b`, `refuge_user`, which is true
for groups utilize that predation refuge, and `bad_pred`, which is false
for predator groups whose body shape or predatory strategy allow them to
access fish within refuge (e.g. eels).

To ensure some food is always available to predators, the maximum
proportion of fish protected by refuge in any size class is set by
`max_protect`.

The refuge profile is used when calculating the food encounter rate in
[`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md)
and the predation mortality rate in
[`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md).
Its entries are dimensionless values between 0 and 1 which represent the
proportion of fish in the corresponding prey and size categories that
are hidden within refuge and thus cannot be encountered by predators. If
no refuge is available then predator-prey interactions are determined
entirely by size-preference.

The mizerReef package provides three methods to define the refuge
profile.

- **Sigmoidal Method**:  
  This method is preferred for data-poor reefs or reefs where the refuge
  distribution is unknown. It is also ideal for systems where only one
  functional group is expected to be utilizing refuge. The proportion of
  fish with access to refuge \\ R_j(w_p) \\ is given by

  \$\$ R_j(w_p) = \frac{r} {1 + e^{ \left( \Delta (w - W\_{refuge})
  \right)} }\$\$

  Here \\W\_{refuge}\\ marks the body weight at which refuge becomes
  scarcer for prey. \\r\\ defines the maximum proportion of fish with
  access to predation refuge and is always less than or equal to
  `max_protect`. \\\alpha\\ controls the rate at which the availability
  of refuge decreases with increasing body size. It defaults to a steep
  slope of 100.

  For this method, `method_params` should contain columns named
  `prop_protect` and `L_refuge` that give the values for \\r\\ and the
  length at which refuge becomes scarce in cm.

- **Binned Method**:  
  This method is appropriate for theoretical applications and does not
  rely on empirical data. It sets refuge to a constant proportion of
  fish within a given size range. The proportion of fish in group \\j\\
  with access to refuge is given by

  \$\$ R_j(w_p) = r_k \~\~\~\~\~\~~ w_p ∈ (~w\_{k-1}, w_k~\] \$\$

  where \\r_k\\ is the proportion of fish with access to refuge in size
  class \\k\\.

  For this method, `method_params` should contain columns named
  `start_L`and `end_L` which contain the starting and ending lengths
  [cm](https://rdrr.io/r/grDevices/cm.html) of each size bin and
  `prop_protect`, the proportion of fish protected within each
  corresponding size bin.

- **Competitive Method**:  
  This method is appropriate when refuge density data is available for
  the modelled reef. The refuge density describes the distribution of
  refuges \\(no./m^2)\\ across defined fish body size categories. The
  proportion of fish in size class \\k\\ with access to refuge is given
  by

  \$\$R\_{j}(w_p) = \tau \cdot \frac{ \eta\_{k} } { \sum_i
  \int\_{w\_{k-1}}^{w-k} N_i(w) \\ dw}\$\$

  where \\ \tau \\ is the proportion of fish with access to refuge that
  are expected to actually utilize it, \\ \eta\_{k}\\ is the density of
  refuges in size range \\(w\_{k-1}, w_k\]\\ and \\\sum\_{i}
  \int\_{w\_{k-1}}^{w_k} N_i(w)~dw\\ gives the density of fish from any
  group in size range \\(w\_{k-1}, w_k\]\\. This represents the density
  of competitors for refuges in size class \\k\\.

  For this method, `method_params` should contain columns named
  `start_L`and `end_L` which contain the starting and ending lengths
  [cm](https://rdrr.io/r/grDevices/cm.html) of each size bin and
  `refuge_density`, the number of refuges available in each size bin
  (no/m^2).

Users can also set a noncomplex reef with no habitat refuge. This option
is convenient for finding steady state parameters.

This function checks that the supplied refuge parameters are valid, adds
relevant columns to the `species_params` data frame, and stores refuge
parameters in the `other_params` slot of the `params` object.

Refuge profile parameters can be input in a spreadsheet program and
saved as a .csv file. The data can then be read into R using the command
[`read.csv()`](https://rdrr.io/r/utils/read.table.html).
