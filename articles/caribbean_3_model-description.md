# Simple 3-Group Example

## Model Description

[Skip to Resources](#resources)

`caribbean_3_model` reproduces the three-group trait-based reef model
from Rogers et al. (2018) – predators, herbivores, and a general group
of benthic invertebrates – inside mizerReef. It is deliberately coarse:
useful when a species-group resolution is enough, or as a minimal
starting point for setting up a new reef model, before working up to
something like the [10-group Karpata Reef
model](https://cmbeese.github.io/mizerReef/articles/karpata_model-description.md),
which parameterises the same three broad functional groups against finer
field data. This page documents this model’s specific parameters and how
they were tuned; for the general mizerReef formulation these parameters
plug into, see
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md).

Biomass estimates were based on data collected from relatively un-fished
reefs within a marine reserve in Bonaire. This data was part of the
FORCE dataset that categorised fish community structure and habitat
complexity across Caribbean coral reefs (Rogers et al. 2014, Williams et
al. 2015, 2016, Newman et al. 2015, Dryden 2016). A marine-reserve site
was chosen deliberately, for this model and for the 10-group Karpata
Reef model alike, so that observed biomasses would not need to be
corrected for fishing pressure before being used to calibrate the
models. Visual surveys recorded fish up to 5 cm in length. The cutoff
weight for observed biomass estimates was calculated using published
length-weight relationships (Froese and Pauly 2023). The cutoff size is
given in the table below along with the observed biomass per square
meter in grams.

| Species Group | Observed Biomass \[g/m^2\] | Cutoff Size \[g\] |
|:--------------|---------------------------:|------------------:|
| predators     |                        107 |             3.125 |
| herbivores    |                         34 |             3.125 |
| inverts       |                         40 |             3.125 |

### Interaction matrix

The $`\theta_{ij}`$ matrix sets the interaction strength between
predator group $`i`$ and prey group $`j`$, with entries between 0
(groups cannot interact) and 1.

![Interaction strength between each predator group (rows) and prey group
(columns). Predators consume all three groups, at different strengths;
herbivores and invertebrates do not prey on
fish.](caribbean_3_model-description_files/figure-html/unnamed-chunk-3-1.png)

Interaction strength between each predator group (rows) and prey group
(columns). Predators consume all three groups, at different strengths;
herbivores and invertebrates do not prey on fish.

### Predation kernel

The parameters for the predation kernels were based on the global
average predator-prey mass ratio for marine organisms and home range
size estimates from Nash et al. (2014). All estimates fall within
observed ranges. All groups use a lognormal predation kernel with the
same parameters as in Rogers et al. (2014).

|            | beta | sigma |     gamma |
|:-----------|-----:|------:|----------:|
| predators  |  100 |     1 | 0.3386893 |
| herbivores |  100 |     1 | 0.4386532 |
| inverts    |  100 |     1 | 1.1214342 |

### Vulnerability to predation

This model uses the refuge method (see
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md)’s
[refuge
profiles](https://cmbeese.github.io/mizerReef/articles/model-description.html#refuge-profiles)
section for the general formulation and Rogers et al. (2018) for its
derivation). All fish are assumed to utilize refuge ($`\tau = 1`$). Fish
smaller than 0.1 g are assumed to be larval reef fish that have not yet
settled to the reef. A maximum proportion of 98 % of fish are protected
at any given time.

Refuge length bins and densities used in the steady state (same Bonaire
FORCE data as the 10-group model)

The following table gives the fish length bins and refuge densities that
define the refuge profile used in the steady state. Length bins and
refuge densities were garnered from the FORCE data in the same location
as biomass estimates – the same
[`karpata_refuge`](https://cmbeese.github.io/mizerReef/reference/karpata_refuge.md)
profile used by the 10-group model.

| Start of Bin (cm) | End of Bin (cm) | Refuge Density (no./m^2) |
|------------------:|----------------:|-------------------------:|
|                 0 |               5 |                    7.533 |
|                 5 |              10 |                    1.400 |
|                10 |              15 |                    0.708 |
|                15 |              20 |                    0.283 |
|                20 |              25 |                    0.100 |
|                25 |              30 |                    0.050 |
|                30 |              35 |                    0.042 |
|                35 |              40 |                    0.033 |
|                40 |              45 |                    0.033 |
|                45 |              50 |                    0.042 |

How each species group interacts with predation refuge

Each species group can utilise benthic structures differently depending
on their specific traits.

| Group      | Uses refuge? | Accesses prey in refuge? |
|:-----------|:-------------|:-------------------------|
| predators  | Yes          | ×                        |
| herbivores | Yes          | Yes                      |
| inverts    | ×            | Yes                      |

The figure below shows the density-dependent refuge profile produced by
the competitive method at steady state. Invertebrates are not shown as
they do not use refuge (`refuge_user = FALSE`, see the table above).

![Proportion of each refuge-using group protected from predation, by
body length, at steady
state.](caribbean_3_model-description_files/figure-html/unnamed-chunk-8-1.png)

Proportion of each refuge-using group protected from predation, by body
length, at steady state.

### Consumption

Consumption of the detrital and algal resources is subject to a Holling
type II functional response (see
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md)’s
[Consumption and
satiation](https://cmbeese.github.io/mizerReef/articles/model-description.html#consumption-and-satiation)
section for the general formula). The parameter $`h`$ is the max
consumption rate for a consumer of size 1 gram, chosen so that consumers
subject to satiation are neither too starved nor totally satiated;
$`\alpha_i`$ is the proportion of consumed biomass retained, with
$`1-\alpha_i`$ expelled as faeces, which contribute to detritus.

No maximum consumption rate is imposed for predators. Herbivores are
subject to satiation here, unlike the package default for herbivorous
groups: recalibrating this model against the corrected
senescence-mortality formula showed herbivore biomass has no
density-dependent brake at all without some cap on individual intake
once mortality is realistically low. The realised feeding level for
herbivores is consistently close to 1, consistent with Caribbean
herbivores’ guts being observed to be full essentially continuously
(Ferreira et al. 1998, Kopp et al. 2010) – the citations behind the
package’s herbivore-satiation default describe grazing *pressure on the
shared algae resource* not self-regulating when food is abundant, which
is unaffected by this change (see [the main model description’s Algae
section](https://cmbeese.github.io/mizerReef/articles/model-description.html#algae)).

|            |        h | alpha |    n |
|:-----------|---------:|------:|-----:|
| herbivores | 11.73461 |   0.6 | 0.75 |
| inverts    | 30.00000 |   0.6 | 0.75 |

#### Metabolic losses

Standard metabolism occurs at rate $`k_{s.i}`$; losses due to activity
and movement occur at rate $`k_i`$ (see
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md)
for the general growth/reproduction formulation).

|            |        ks |    p |   k |
|:-----------|----------:|-----:|----:|
| predators  | 0.0543037 | 0.75 |   0 |
| herbivores | 0.0468876 | 0.75 |   0 |
| inverts    | 0.1000000 | 0.75 |   0 |

#### Maturation and growth

Maturation length and age data were based on the most commonly observed
species from each functional group in the FORCE data set (Rogers et al.
2014, Williams et al. 2015, 2016, Newman et al. 2015, Dryden 2016). For
predators, this was the graysby grouper, *Cephalopholis cruentata*. The
most commonly observed herbivore was the stoplight parrotfish,
*Sparisoma viride*. Relevant maturity parameters were pulled from
FishBase (Froese and Pauly 2023).

The herbivore `age_mat` value is set to 1.6 years, the age at median
sexual maturity (AM50) reported for *S. viride* by Hernández and
Shervette (2025) – a 2013-2023 otolith/gonad-histology study of 1801
U.S. Caribbean stoplight parrotfish. An earlier version of this model
used `age_mat = 4`, which appears to have conflated AM50 with that same
study’s age at median sexual transition (AT50 = 4.5 years, the age at
which female stoplight parrotfish – a protogynous hermaphrodite species
– become male), a distinct milestone from first reproductive maturity.

|            | w_mat | w_max | age_mat |
|:-----------|------:|------:|--------:|
| predators  | 102.4 |  3125 |       2 |
| herbivores | 102.4 |  3125 |       2 |
| inverts    |   0.1 |  3125 |      NA |

The values for the growth parameters below were chosen so that the
resulting growth curves would be close to von Bertalanffy growth curves;
$`a`$ and $`b`$ are the allometric weight-length parameters
$`w = a l^b`$ ($`w`$ in grams, $`l`$ in centimetres), taken from the
literature.

|            | k_vb | w_max |     a |   b |
|:-----------|-----:|------:|------:|----:|
| predators  |  0.4 |  3125 | 0.025 |   3 |
| herbivores |  0.6 |  3125 | 0.025 |   3 |
| inverts    |   NA |  3125 | 0.025 |   3 |

### Mortality

See
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md)’s
[Mortality](https://cmbeese.github.io/mizerReef/articles/model-description.html#mortality)
section for the general formulation. This model’s residual-natural- and
senescence-mortality parameters:

We use a residual natural mortality rate of $`\mu_{nat} =`$ 0.2 per year
at size 1 gram, and senescence mortality parameters $`k_{sen} =`$ 0.1
(`sen_prop`) and $`p_{sen} =`$ 0.3 (`sen_curve`, the exponent governing
how steeply mortality climbs as individuals approach their maximum
size), both based on estimates from Hatcher (1988).

#### Fishing mortality

Fishing parameters for mizerReef models can be set up with
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.html).
The vignettes were run with the following fishing parameters:

| Group      | Minimum fishing size \[g\] | Catchability \[1/year\] |
|:-----------|---------------------------:|------------------------:|
| predators  |                      102.4 |                       1 |
| herbivores |                      102.4 |                       1 |

### Reproduction

The reproduction parameters $`\epsilon_i`$ and $`R_{max.i}`$ are not
directly observable. Their values were instead chosen so as to produce
steady-state abundances of the groups that are in line with observations
and to give reasonable values for the reproduction level – the ratio
between the actual reproduction rate $`R_i`$ and the maximal possible
reproduction rate $`R_{\max.i}`$.

|            | w_min |   erepro |      R_max |
|:-----------|------:|---------:|-----------:|
| predators  | 0.001 | 2.03e-05 |  0.1662314 |
| herbivores | 0.001 | 4.00e-07 |  0.0263658 |
| inverts    | 0.001 | 8.28e-04 | 93.0679393 |

## Resources

[Skip to Tuning the Steady State](#tuning-the-steady-state)

The fish spectrum is fed by three background resources: plankton, algae,
and detritus. Small individuals of all species feed on plankton, while
only herbivorous groups and invertebrates feed on algae or detritus. See
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md)’s
[unstructured resource
dynamics](https://cmbeese.github.io/mizerReef/articles/model-description.html#unstructured-resource-dynamics)
section for the general formulation and its relationship to mizerShelf.

#### Plankton

The plankton spectrum ranges from $`w_0=9\times 10^{-13}`$ to
$`w_{cutoff}=1`$ grams. The steady state abundance of plankton at size 1
gram is $`\kappa=11.8 (g/m^{2})`$, with slope $`\lambda=2.05`$.

#### Algae

The algal resource is described only by its total biomass $`B_A`$, and
feeding on it is not size-based. In the steady state the total algal
biomass per square meter is $`B_A = 350`$ grams.

|            |       rho | interaction_algae |
|:-----------|----------:|------------------:|
| herbivores | 0.4840935 |                 1 |

Algal *production* – a fixed, literature-informed constant,
`algae_growth_initial` in
[`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)
– is 2000 grams per square meter per year at steady state. Note this
value is not retuned to match consumption the way detritus production is
(see below); see
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md)’s
Algae section for why.

#### Detritus

Detritus is consumed by herbivores and benthic invertebrates, also not
size-based. In the steady state the total detrital biomass per square
meter is $`B_D = 0.009575`$ grams.

|         |      rho | interaction_detritus |
|:--------|---------:|---------------------:|
| inverts | 276.4376 |                    1 |

Detritus production comes from three sources: defecation, decomposing
dead organisms, and external input. At steady state, decomposing dead
organisms contribute 15.05 g/year and defecation contributes 56.66
g/year – 20 % of external mortality and 80 % of senescence mortality is
assumed to decompose to detritus (estimates from Hatcher (1988)). The
remaining external-input term is solved for so that, unlike algae,
detritus *production* matches current consumption: at steady state it is
-49.93 g/year. Where this is negative, faeces and mortality alone
already produce more detritus than is consumed, and the excess is
assumed to be washed away.

## Tuning the Steady State

The following R script was used to tune the steady state parameters for
this model (see `inst/scripts/Caribbean_3_model-calibration.R` in the
package source for the exact, currently-runnable version, and
[`vignette("steady-state-recipe")`](https://cmbeese.github.io/mizerReef/articles/steady-state-recipe.md)
for the general recipe this follows, including the final algae/detritus
absolute-scale rescaling step).

Show the full steady-state calibration script

``` r

# Setting up a generic Caribbean coral reef model with multiple resources
# Three groups: Predators, Herbivores, Invertebrates
# Model steady state calibration

## Setup - load packages --------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(here)

## Load parameters ----------------------------------------------------------
caribbean_3_species <- read.csv(here("inst/data-csv/caribbean_3_species.csv"))
caribbean_3_interaction <- read.csv(here("inst/data-csv/caribbean_3_interaction.csv"),
                                    row.names = 1)
# Refuge densities from Karpata reef in Bonaire, FORCE dataset
karpata_refuge <- read.csv(here("inst/data-csv/karpata_refuge.csv"))
tuning_profile <- read.csv(here("inst/data-csv/tuning_profile.csv"))

# With these parameters, herbivores consume plankton at small sizes and
#   transition fully to algae by maturity
# With these parameters, invertebrates consume plankton and detritus,
#   with the proportion of detritus increasing with size

## Set model, with info_level = 1 to keep only the reports that say
## something went differently from what was asked, not every default fill ---
params <- newReefParams(
    species_params = caribbean_3_species,
    interaction = caribbean_3_interaction,
    method = "binned",
    method_params = tuning_profile,
    info_level = 1
)

## Project to first steady state ----------------------------------------------
params <- reefSteady(params)

## Calibrate biomasses and growth ----------------------------------------------

# Match observed species group biomasses
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- reefSteady(params)

# Match observed growth rates
params <- matchReefGrowth(params)
params <- reefSteady(params)

# Iterate to refine biomass - repeat many times; piping keeps it readable.
params <- params |>
    calibrateReefBiomass() |> matchBiomasses() |> matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses() |> matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses() |> matchReefGrowth() |>
    reefSteady()
# (repeated several more times in the full calibration script until
#  plotBiomassVsSpecies(params) matches observations closely)

plotBiomassVsSpecies(params) # spot on
plotTotalAbundance(params)
plotTotalBiomass(params)

## Now switch to competitive method ---------------------------------------------
params <- newRefuge(params,
    new_method = "competitive",
    new_method_params = karpata_refuge
)

# Match biomasses again with the new refuge profile
params <- params |>
    calibrateReefBiomass() |> matchBiomasses() |> matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |> matchBiomasses() |> matchReefGrowth() |>
    reefSteady()
# (repeated several more times in the full calibration script)

# Make sure new refuge is in place
plotVulnerable(params)
plotRefugeProfile(params)
plotBiomassVsSpecies(params) # spot on

## Check resulting spectra and tune resources ------------------------------------

# Spectra should be reasonably straight to match predictions of Sheldon's
# spectrum but also have nonlinearities at refuge sizes
plotSpectra(params, total = TRUE, biomass = TRUE)
plotSpectra(params, total = TRUE, biomass = TRUE, per_log_size = TRUE)
plotFeedingLevel(params, species = "inverts")

# Tune reproduction -----------------------------------------------------------
# We do not have yield or catch data - can't tune size distribution
params <- setBevertonHolt(params, erepro = 0.0134)
params <- reefSteady(params)

# Increase reproduction level to reduce excessive density dependence
rep_level <- c(0.5, 0.5, getReproductionLevel(params)["inverts"])
names(rep_level) <- c("predators", "herbivores", "inverts")
params <- setBevertonHolt(params, reproduction_level = rep_level)

# Iterate to get back to steady state
params <- params |> reefSteady() |> reefSteady() |> reefSteady()

# Rescale algae/detritus to a realistic absolute scale, now that diet,
# biomass and growth are tuned -- see vignette("steady-state-recipe")'s
# final step for the literature targets and why this belongs last.
params <- rescale_algae(params, target_biomass / algae_biomass(params))
detritus_lifetime(params) <- target_lifetime

# Final diagnostic plots
plotBiomassVsSpecies(params)
plotRefugeProfile(params)
plotSpectra(params, biomass = TRUE, total = TRUE)
plotDiet(params)
plotGrowthCurves(params)
plotPredMort(params)

# Save!
caribbean_3_model <- reefSteady(params)
save(caribbean_3_model, file = "data/caribbean_3_model.rda")
```

Dryden, C. 2016. Habitat structural complexity of caribbean coral reefs
and its relationships with fish community structure. PhD thesis,
Newcastle University, Newcastle upon Tyne United Kingdom.

Ferreira, D. E. L., A. C. Peret, and R. Coutinho. 1998. [Seasonal
grazing rates and food processing by tropical herbivorous
fishes](https://doi.org/10.1111/j.1095-8649.1998.tb01029.x). Journal of
Fish Biology 53:222–235.

Froese, R., and D. Pauly. 2023. [FishBase](https://www.fishbase.org).

Hatcher, B. G. 1988. [Coral reef primary productivity: A beggar’s
banquet](https://doi.org/10.1016/0169-5347(88)90117-6). Trends in
Ecology & Evolution 3:106–111.

Hernández, J. M. R., and V. R. Shervette. 2025. [Addressing life history
information gaps for Caribbean parrotfishes: Queen parrotfish Scarus
vetula and stoplight parrotfish Sparisoma
viride](https://doi.org/10.1007/s10641-024-01651-x). Environmental
Biology of Fishes 108:179–198.

Kopp, D., Y. Bouchon-Navaro, S. Cordonnier, A. Haouisée, M. Louis, and
C. Bouchon. 2010. [Evaluation of algal regulation by herbivorous fishes
on caribbean coral reefs](https://doi.org/10.1007/s10152-009-0177-4).
Helgoland Marine Research 64:181–190.

Nash, K. L., J. Q. Welsh, N. A. J. Graham, and D. R. Bellwood. 2014.
[Home-range allometry in coral reef fishes: Comparison to other
vertebrates, methodological issues and management
implications](https://doi.org/10.1007/s00442-014-3152-y). Oecologia
177:73–83.

Newman, S. P., E. H. Meesters, C. S. Dryden, S. M. Williams, C. Sanchez,
P. J. Mumby, and N. V. C. Polunin. 2015. [Reef flattening effects on
total richness and species responses in the
caribbean](https://doi.org/10.1111/1365-2656.12429). Journal of Animal
Ecology 84:1678–1689.

Rogers, A., J. L. Blanchard, and P. J. Mumby. 2014. [Vulnerability of
coral reef fisheries to a loss of structural
complexity](https://doi.org/10.1016/j.cub.2014.03.026). Current Biology
24:1000–1005.

Rogers, A., J. L. Blanchard, S. P. Newman, C. S. Dryden, and P. J.
Mumby. 2018. [High refuge availability on coral reefs increases the
vulnerability of reef‐associated predators to
overexploitation](https://doi.org/10.1002/ecy.2103). Ecology 99:450–463.

Williams, S. M., I. Chollett, G. Roff, J. Cortaes, C. S. Dryden, and P.
J. Mumby. 2015. [Hierarchical spatial patterns in caribbean reef benthic
assemblages](https://doi.org/10.1111/jbi.12509). Journal of Biogeography
42:1327–1335.

Williams, S. M., C. S’anchez?God’ınez, S. P. Newman, and J. Cort’es.
2016. [Ecological assessments of the coral reef communities in the
eastern caribbean and the effects of herbivory in influencing coral
juvenile density and algal cover](https://doi.org/10.1111/maec.12395).
Marine Ecology 38.
