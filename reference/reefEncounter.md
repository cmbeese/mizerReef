# Get encounter rate needed to project a mizerReef model

Calculates the rate \\E_i(w)\\ at which a predator from group \\i\\ and
weight \\w\\ encounters food (grams/year). You would not usually call
this function directly but instead use
[`getEncounter()`](https://sizespectrum.org/mizer/reference/getEncounter.html),
which then calls this function.

## Usage

``` r
reefEncounter(params, n, n_pp, n_other, t, vulnerable = NULL, ...)
```

## Arguments

- params:

  A
  [MizerParams](https://sizespectrum.org/mizer/reference/MizerParams.html)
  object

- n:

  A matrix of species abundances (species x size).

- n_pp:

  A vector of the resource abundance by size

- n_other:

  A list of abundances for other dynamical components of the ecosystem

- t:

  The time for which to do the calculation (Not used by standard mizer
  rate functions but useful for extensions with time-dependent
  parameters.)

- vulnerable:

  A two dimensional array (prey species x prey size) with the proportion
  of prey vulnerable to being encountered.

- ...:

  Unused

## Value

A named two dimensional array (predator species x predator size) with
the encounter rates.

## Predation encounter

The encounter rate \\E_i(w)\\ at which a predator of species \\i\\ and
weight \\w\\ encounters food has contributions from the encounter of
fish prey and of resources. This is determined by summing over all prey
species and the resource spectrum and then integrating over all prey
sizes \\w_p\\, weighted by predation kernel \\\phi(w,w_p)\\:

\$\$ E_i(w) = \gamma_i(w) \int \left( \theta\_{ip} N_R(w_p) + \sum\_{j}
\theta\_{ij} V\_{ji}(w_p) N_j(w_p) \right) \phi_i(w,w_p) w_p \\ dw_p.
\$\$ Here \\N_j(w)\\ is the abundance density of species \\j\\ and
\\N_R(w)\\ is the abundance density of resource. The overall prefactor
\\\gamma_i(w)\\ determines the predation power of the predator. It could
be interpreted as a search volume and is set with the
[`setSearchVolume()`](https://sizespectrum.org/mizer/reference/setSearchVolume.html)
function.

The predation kernel \\\phi(w,w_p)\\is set with the
[`setPredKernel()`](https://sizespectrum.org/mizer/reference/setPredKernel.html)
function.

The vulnerability to predation, \\V\_{ji}(w)\\ accounts for protective
behavior of the prey. The parameters that control this are set with the
[`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)
function.

The species interaction matrix \\\theta\_{ij}\\ is set with
[`setInteraction()`](https://sizespectrum.org/mizer/reference/setInteraction.html)
and the resource interaction vector \\\theta\_{ip}\\ is taken from the
`interaction_resource`column in `params@species_params`.

## Details

The encounter rate is multiplied by \\1-f_0\\ to obtain the consumption
rate, where \\f_0\\ is the feeding level calculated with
[`getFeedingLevel()`](https://sizespectrum.org/mizer/reference/getFeedingLevel.html).
This is used by the
[`project()`](https://sizespectrum.org/mizer/reference/project.html)
function for performing simulations.

The function returns values also for sizes outside the size-range of the
species. These values should not be used, as they are meaningless.

If your model contains additional components that you added with
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)
and for which you specified an `encounter_fun` function then the
encounters of these components will be included in the returned value.

## See also

Other mizer rate functions:
[`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md),
[`reefFeedingLevel()`](https://cmbeese.github.io/mizerReef/reference/reefFeedingLevel.md),
[`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md),
[`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md),
[`reefRates()`](https://cmbeese.github.io/mizerReef/reference/reefRates.md),
[`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
