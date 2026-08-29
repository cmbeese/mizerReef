# Get vulnerability level at in time range t

Returns the proportion of fish at size \\w\\ that are not hidden in
predation refuge and thus vulnerable to being encountered by predators.

## Usage

``` r
getDegrade(object, n, n_pp, n_other, time_range, drop = TRUE, ...)
```

## Arguments

- object:

  A `MizerParams` object or a `MizerSim` object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components of the ecosystem

- time_range:

  A numeric or character vector of times. Only the range of values
  matters, so all saved times between `min(time_range)` and
  `max(time_range)` are selected.

- drop:

  If `TRUE` then any dimension of length 1 will be removed from the
  returned array.

- ...:

  Unused

## Value

If a `MizerParams` object is passed in, the function returns a numeric
vector of the refuge density for each size bin.

         If a `MizerSim` object is passed in, the function returns a two
         dimensional array (time step x refuge size bin) with the refuge
         density calculated at every time step in the simulation. If
         \code{drop = TRUE} then the dimension of length 1 will be
         removed from the returned array.

## Details

This function uses
[`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
to calculate the vulnerability to predation.

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

Other rate functions:
[`getSenMort()`](https://cmbeese.github.io/mizerReef/reference/getSenMort.md),
[`getVulnerable()`](https://cmbeese.github.io/mizerReef/reference/getVulnerable.md)
