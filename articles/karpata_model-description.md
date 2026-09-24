# Functional Group Example

## Model Description

[Skip to Resources](#resources)

`caribbean_10_model` is mizerReef’s field-calibrated flagship example: 9
fish species groups plus a general group of benthic invertebrates,
parameterised against field data from Karpata Reef, Bonaire. Species
groups were assigned based on functional traits, body size, diet, and
interactions with habitat structure. It is the finer-resolution
counterpart to the [3-group Caribbean
model](https://cmbeese.github.io/mizerReef/articles/caribbean_3_model-description.md),
which reproduces the coarser trait-based structure of Rogers et al.
(2018) using the same marine-reserve dataset described below. This page
documents this model’s specific parameters and how they were tuned; for
the general mizerReef formulation these parameters plug into, see
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md).

Biomass estimates were based on data collected from Karpata Reef in
Bonaire, a site with low fishing levels in a marine reserve (Rogers et
al. 2014, Williams et al. 2015, 2016, Newman et al. 2015, Dryden 2016).
Observations of fish smaller than 10 cm were removed from analysis
because their abundance is not well captured by visual survey methods.
The cutoff weight for observed biomass estimates was calculated using
published length-weight relationships (Froese and Pauly 2023). The
cutoff size is given in the table below along with the observed biomass
per square meter in grams.

| Species Group | Observed Biomass \[g/m^2\] | Cutoff Size \[g\] |
|:--------------|---------------------------:|------------------:|
| pred_eng      |                       5.51 |          12.62969 |
| pred_grab     |                      27.60 |          17.80530 |
| eels          |                      10.00 |                NA |
| pred_inv      |                      14.13 |          26.58539 |
| pred_plank    |                       0.55 |          13.49043 |
| parrotfish    |                      30.56 |          15.48385 |
| farm_damsel   |                       0.40 |                NA |
| herbs         |                       1.50 |          75.75344 |

### Interaction matrix

The $`\theta_{ij}`$ matrix sets the interaction strength between
predator group $`i`$ and prey group $`j`$, with entries between 0
(groups cannot interact) and 1. All organisms were assumed to interact
equally, with the exception of nocturnal invertivores which are known to
hunt off the reef at night.

![Interaction strength between each predator group (rows) and prey group
(columns) in the 10-group Karpata Reef
model.](karpata_model-description_files/figure-html/unnamed-chunk-3-1.png)

Interaction strength between each predator group (rows) and prey group
(columns) in the 10-group Karpata Reef model.

### Predation kernel

The parameters for the predation kernels were based on trait-based
studies of prey size selectivity, diet studies that estimate
predator-prey mass ratio, and home range size estimates from Nash et al.
(2014), tuned to achieve expected diet compositions. All estimates fall
within observed ranges. All groups use a lognormal predation kernel.

|             | beta | sigma |     gamma |
|:------------|-----:|------:|----------:|
| pred_eng    |   50 |     1 | 1.9889276 |
| pred_grab   |   30 |     1 | 1.6527159 |
| eels        |   30 |     2 | 0.2180816 |
| pred_crypt  |    5 |     1 | 2.5541765 |
| pred_inv    |   30 |     2 | 1.9941828 |
| pred_plank  | 5000 |     2 | 0.2102966 |
| parrotfish  |   30 |     1 | 2.5999971 |
| farm_damsel |   30 |     1 | 0.6117167 |
| herbs       |   30 |     1 | 2.4050383 |
| inverts     |   30 |     1 | 9.6181106 |

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

Refuge length bins and densities used in the steady state (Karpata Reef,
FORCE data)

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

Each species group utilises benthic structures differently depending on
their specific traits.

| Group       | Uses refuge? | Accesses prey in refuge? |
|:------------|:-------------|:-------------------------|
| pred_eng    | Yes          | ×                        |
| pred_grab   | Yes          | ×                        |
| eels        | Yes          | Yes                      |
| pred_crypt  | Yes          | Yes                      |
| pred_inv    | Yes          | Yes                      |
| pred_plank  | Yes          | Yes                      |
| parrotfish  | Yes          | Yes                      |
| farm_damsel | Yes          | Yes                      |
| herbs       | Yes          | Yes                      |
| inverts     | ×            | Yes                      |

The figure below shows the density-dependent refuge profile produced by
the competitive method at steady state.

![Proportion of each refuge-using group protected from predation, by
body length, at steady state, for the 10-group Karpata Reef
model.](karpata_model-description_files/figure-html/unnamed-chunk-8-1.png)

Proportion of each refuge-using group protected from predation, by body
length, at steady state, for the 10-group Karpata Reef model.

### Consumption

Invertebrate consumption of the detrital resource and planktivory are
subject to a Holling type II functional response (see
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md)’s
[Consumption and
satiation](https://cmbeese.github.io/mizerReef/articles/model-description.html#consumption-and-satiation)
section for the general formula). $`h`$ is the max consumption rate for
an invertebrate consumer of size 1 gram, chosen so that invertebrates
are neither too starved nor totally satiated; $`\alpha_i`$ is the
proportion of consumed biomass retained, with $`1-\alpha_i`$ expelled as
faeces, which contribute to detritus. No maximum consumption rate is
imposed for predatory or herbivorous groups.

|             |         h | alpha |    n |
|:------------|----------:|------:|-----:|
| pred_plank  |  13.36972 |   0.6 | 0.75 |
| parrotfish  |  63.73181 |   0.6 | 0.75 |
| farm_damsel |  14.99456 |   0.6 | 0.75 |
| herbs       |  58.95293 |   0.6 | 0.75 |
| inverts     | 235.76163 |   0.6 | 0.75 |

#### Metabolic losses

Standard metabolism occurs at rate $`k_{s.i}`$; losses due to activity
and movement occur at rate $`k_i`$ (see
[`vignette("model-description")`](https://cmbeese.github.io/mizerReef/articles/model-description.md)
for the general growth/reproduction formulation).

|             |        ks |    p |   k |
|:------------|----------:|-----:|----:|
| pred_eng    | 0.3629290 | 0.75 |   0 |
| pred_grab   | 0.2235502 | 0.75 |   0 |
| eels        | 0.0699183 | 0.75 |   0 |
| pred_crypt  | 0.1500000 | 0.75 |   0 |
| pred_inv    | 0.1710873 | 0.75 |   0 |
| pred_plank  | 0.1114551 | 0.75 |   0 |
| parrotfish  | 0.2897465 | 0.75 |   0 |
| farm_damsel | 0.1094264 | 0.75 |   0 |
| herbs       | 0.2339976 | 0.75 |   0 |
| inverts     | 0.1000000 | 0.75 |   0 |

#### Maturation and growth

Maturation length and age data were based on the most observed species
from each species group in the FORCE data set (Rogers et al. 2014,
Williams et al. 2015, 2016, Newman et al. 2015, Dryden 2016).

|             | w_mat |       w_max | age_mat |
|:------------|------:|------------:|--------:|
| pred_eng    |  50.0 | 1259.571941 |     2.0 |
| pred_grab   | 140.0 | 3915.475510 |     5.4 |
| eels        | 300.0 | 2959.552686 |     3.0 |
| pred_crypt  |   0.8 |    6.242901 |      NA |
| pred_inv    |  50.0 | 1600.073891 |     0.5 |
| pred_plank  |   2.5 |  110.191125 |     1.0 |
| parrotfish  |  63.0 | 4453.772271 |     1.6 |
| farm_damsel |   1.0 |   41.540096 |     1.0 |
| herbs       | 105.0 | 1269.327573 |     2.0 |
| inverts     |   0.1 |  675.000000 |      NA |

The values for the growth parameters below were chosen so that the
resulting growth curves would be close to von Bertalanffy growth curves;
$`a`$ and $`b`$ are the allometric weight-length parameters
$`w = a l^b`$ ($`w`$ in grams, $`l`$ in centimetres), taken from the
literature. There is no negative growth or starvation mortality.

|             | k_vb |       w_max |       a |    b |
|:------------|-----:|------------:|--------:|-----:|
| pred_eng    | 0.40 | 1259.571941 | 0.01100 | 3.06 |
| pred_grab   | 0.10 | 3915.475510 | 0.01740 | 3.01 |
| eels        | 0.19 | 2959.552686 | 0.00098 | 3.24 |
| pred_crypt  | 3.00 |    6.242901 | 0.01122 | 3.04 |
| pred_inv    | 0.40 | 1600.073891 | 0.01200 | 3.10 |
| pred_plank  | 0.30 |  110.191125 | 0.01259 | 3.03 |
| parrotfish  | 0.60 | 4453.772271 | 0.01380 | 3.05 |
| farm_damsel | 0.30 |   41.540096 | 0.02042 | 2.97 |
| herbs       | 0.40 | 1269.327573 | 0.02570 | 2.95 |
| inverts     | 2.00 |  675.000000 | 0.02500 | 3.00 |

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

| Group       | Minimum fishing size \[g\] | Catchability \[1/year\] |
|:------------|---------------------------:|------------------------:|
| pred_eng    |                       50.0 |                       1 |
| pred_grab   |                      140.0 |                       1 |
| eels        |                      300.0 |                       1 |
| pred_crypt  |                        0.8 |                       1 |
| pred_inv    |                       50.0 |                       1 |
| pred_plank  |                        2.5 |                       1 |
| parrotfish  |                       63.0 |                       1 |
| farm_damsel |                        1.0 |                       1 |
| herbs       |                      105.0 |                       1 |
| inverts     |                        0.1 |                       1 |

### Reproduction

The reproduction parameters $`\epsilon_i`$ and $`R_{max.i}`$ are not
directly observable. Their values were instead chosen so as to produce
steady-state abundances of the groups that are in line with observations
and to give reasonable values for the reproduction level – the ratio
between the actual reproduction rate $`R_i`$ and the maximal possible
reproduction rate $`R_{\max.i}`$.

|             | w_min |     erepro |      R_max |
|:------------|------:|-----------:|-----------:|
| pred_eng    | 0.001 |  0.0078512 |  0.9775429 |
| pred_grab   | 0.001 |  0.0127698 |  9.6063695 |
| eels        | 0.001 |  0.0105685 |  1.9333590 |
| pred_crypt  | 0.001 |  0.0030568 | 36.1032036 |
| pred_inv    | 0.001 |  0.0002371 |  0.5990943 |
| pred_plank  | 0.001 |  0.0751062 | 20.9574142 |
| parrotfish  | 0.001 |  0.0014714 | 15.3921490 |
| farm_damsel | 0.001 | 12.6637873 | 24.7707790 |
| herbs       | 0.001 |  0.0021280 |  0.9593467 |
| inverts     | 0.001 |  0.0090392 |  0.3665275 |

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
$`w_{cutoff}=0.1`$ grams. The steady state abundance of plankton at size
1 gram is $`\kappa=11.4 (g/m^{2})`$, with slope $`\lambda=2.05`$.

#### Algae

The algal resource is described only by its total biomass $`B_A`$, and
feeding on it is not size-based. In the steady state the total algal
biomass per square meter is $`B_A = 4.15\times 10^{-9}`$ grams.

|             |         rho | interaction_algae |
|:------------|------------:|------------------:|
| parrotfish  | 40649043272 |               0.5 |
| farm_damsel |  9563740389 |               0.5 |
| herbs       | 37601006284 |               0.5 |

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
meter is $`B_D = 1.143\times 10^{-10}`$ grams.

|             |          rho | interaction_detritus |
|:------------|-------------:|---------------------:|
| pred_crypt  |  34409747249 |                  0.5 |
| parrotfish  |  40649043272 |                  0.5 |
| farm_damsel |   9563740389 |                  0.5 |
| herbs       |  37601006284 |                  0.5 |
| inverts     | 300744178552 |                  1.0 |

Detritus production comes from three sources: defecation, decomposing
dead organisms, and external input. At steady state, decomposing dead
organisms contribute 46.09 g/year and defecation contributes 284.23
g/year – 80 % of external mortality and 80 % of senescence mortality is
assumed to decompose to detritus (estimates from Hatcher (1988)). The
remaining external-input term is solved for so that, unlike algae,
detritus *production* matches current consumption: at steady state it is
-305.69 g/year. Where this is negative, faeces and mortality alone
already produce more detritus than is consumed, and the excess is
assumed to be washed away.

## Tuning the Steady State

The following R script was used to tune the steady state parameters for
this model (see
[`vignette("steady-state-recipe")`](https://cmbeese.github.io/mizerReef/articles/steady-state-recipe.md)
for the general recipe this follows, including the final algae/detritus
absolute-scale rescaling step).

Show the full steady-state calibration script

``` r

# Setting up a Caribbean coral reef mizer model with multiple resources
# Model steady state calibration
## Setup - load packages ----------------------------------
library(ggplot2)
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(assertthat)
library(here)

## Load parameters -----------------------------------
karpata_10plus  <- read.csv(here("inst/data-csv/caribbean_10_species.csv"))
karpata_int     <- read.csv(here("inst/data-csv/caribbean_10_interaction.csv"),
                            row.names = 1)
karpata_refuge  <- karpata_refuge
tuning_profile  <- tuning_profile

# Herbivores consume plankton at small sizes and
#   transition to detritus and algae as they grow
# Invertebrates consume plankton and detritus,
#   with the proportion of detritus increasing with size

## Set model, with info_level = 1 to keep only the reports that say
## something went differently from what was asked, not every default fill ---
params <- newReefParams(species_params = karpata_10plus,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile,
                        info_level = 1)

## Reduce density dependent of reproduction ----------------
rdi <- rep(0.5, dim(karpata_int)[1])

params <- setBevertonHolt(params, reproduction_level = rdi)
getReproductionLevel(params)

## Project to first steady state -------------------------------
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady()

## Calibrate biomasses and growth ---------------------------------
# Match observed species group biomasses
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- reefSteady(params)

# Match observed growth rates
params <- matchReefGrowth(params)
params <- reefSteady(params)

# Iterate to refine biomass
params <- params |>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|>
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|>
    reefSteady()|>
    calibrateReefBiomass() |> matchBiomasses()|> matchReefGrowth()|>
    reefSteady()

# Check biomass match
plotBiomassVsSpecies(params) # spot on

# Check match with observed age at maturity
age_mat_observed = karpata_10plus$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Closer than needed

# Check predation mortality, feeding levels, and diets
plotPredMort(params) + facet_wrap(~Species)
plotFeedingLevel(params)
plotDiet(params) + scale_x_log10(limits = c(1, 1e4))
plotSpectra(params, biomass = TRUE)

## Now switch to competitive method ----------------------------
params <- newRefuge(params,
                    new_method = "competitive",
                    new_method_params = karpata_refuge)

# Match biomasses again
params <- params |>
    matchBiomasses()|> reefSteady()|>
    matchBiomasses()|> reefSteady()|>
    matchBiomasses()|> reefSteady()|>
    matchBiomasses()|> reefSteady()

# Make sure new refuge is in place
plotVulnerable(params)

plotBiomassVsSpecies(params) # spot on

# Check match with observed age at maturity
age_mat_observed = karpata_10plus$age_mat
age_mat_model = age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Still look good

## Check resulting spectra and tune resources------------------------

# Resource looks low - should match sheldon's spectrum
# looks fairly straight not bad but some bumps
plotSpectra(params, total = TRUE, biomass = TRUE)
plotSpectra(params, total = TRUE, biomass = TRUE, per_log_size = TRUE)

# plot feeding level to check if resource is too low
plotFeedingLevel(params, species = "inverts")

# Invert feeding level is relatively stable through life, non-linearities
#   are probably due to refuge

# Tune reproduction -----------------------------------------------
# We do not have yield or catch data - can't tune size distribution
# First attempt to set very low to see what the minimum values are
params <- setBevertonHolt(params, erepro = 0.0001)
# Now set setting erepro same for all species, as low as possible
params <- setBevertonHolt(params, erepro = 0.35)
# Project back to steady
params <- reefSteady(params)
# Check reproduction level (value between 0 and 1) - should be higher for
# larger, slow growing species and low for small, fast growing ones
getReproductionLevel(params)
# A reproduction level closer to one means reproduction rate is
# almost totally independent of the investment into reproduction
# These are near one for all species except farming damsels

# Check comparison of density dependent & independent reduction
getRDI(params) / getRDD(params)
# Reproduction is density independent for almost all species

# Let's increase reproduction level back to 0.5 so there is still some
# density dependence
params <- setBevertonHolt(params, reproduction_level = rdi)

# Iterate to get back to steady state
params <- params |>
    reefSteady()|>
    reefSteady()|>
    reefSteady()

# Check new reproduction - these look better
rep <- getReproductionLevel(params)
getRDI(params) / getRDD(params)

# Check new spectra & plots
plotSpectra(params, total = TRUE, biomass = TRUE, per_log_size = TRUE)
plotPredMort(params) + facet_wrap(~Species)
plotFeedingLevel(params)
plotDiet(params) + scale_x_log10(limits = c(1, 1e4))
plotSpectra(params, biomass = TRUE)

# Rescale algae/detritus to a realistic absolute scale, now that diet,
# biomass and growth are tuned -- see vignette("steady-state-recipe")'s
# final step for the literature targets and why this belongs last.
params <- rescale_algae(params, target_biomass / algae_biomass(params))
detritus_lifetime(params) <- target_lifetime

# Save!
caribbean_10_model <- reefSteady(params)
```

Dryden, C. 2016. Habitat structural complexity of caribbean coral reefs
and its relationships with fish community structure. PhD thesis,
Newcastle University, Newcastle upon Tyne United Kingdom.

Froese, R., and D. Pauly. 2023. [FishBase](https://www.fishbase.org).

Hatcher, B. G. 1988. [Coral reef primary productivity: A beggar’s
banquet](https://doi.org/10.1016/0169-5347(88)90117-6). Trends in
Ecology & Evolution 3:106–111.

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
