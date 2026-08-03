# Functional Group Example

## Model Description

[Skip to Size Spectrum Dynamics](#size-spectrum-dynamics)

The model includes 9 species groups of fish as well as a general group
of benthic invertebrates. Species groups were assigned based on
functional traits, body size, diet, and interactions with habitat
structure. Estimates of observed biomass for each group were used to
tune reproduction and resource consumption parameters so that steady
state abundances agree with empirical observations. The process of
establishing a steady state that agrees with empirical observations is
nontrivial. The R script used to tune the steady state parameters is
included at the end of this file.

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

## Size Spectrum Dynamics

[Skip to Resources](#resources)

These dynamics expand on those used in standard Mizer models, see Delius
et al. (2023).

### Growth

The growth of organisms is dependent on the energy they are able to
obtain from consumed food resources.

#### Predator-prey encounter rate

The rate $`E_{i}(w)`$ at which a predator of species $`i`$ and weight
$`w`$ encounters food is dictated by the predation power
$`\gamma_i(w)`$, interaction matrix $`\theta_{ij}`$, the vulnerability
to predation $`V_{ij}(w_p)`$, and the size selectivity of predator,
given by the predation kernel $`\phi_i(w,w_p)`$.

##### Interaction Matrix

The $`\theta_{ij}`$ matrix sets the interaction strength between
predator group $`i`$ prey group $`j`$. The predator/prey interaction
matrix has entries between 0 (if the groups can not interact) and 1, see
the figure below. All organisms were assumed to interact equally, with
the exception of nocturnal invertivores which are known to hunt off the
reef at night.

![Interaction strength between each predator group (rows) and prey group
(columns) in the 10-group Karpata Reef
model.](karpata_model-description_files/figure-html/unnamed-chunk-3-1.png)

Interaction strength between each predator group (rows) and prey group
(columns) in the 10-group Karpata Reef model.

##### Predation Kernel

The parameters for the predation kernels were based on trait-based
studies of prey size selectivity, diet studies that estimate predator
prey mass ratio, and home range size estimates from Nash et al. (2014).
Values were tuned to achieve expected diet compositions. All estimates
fall within observed ranges. All groups use a lognormal predation
kernel. The parameters are given in the table below.

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

#### Vulnerability to Predation

The refuge function, $`R_j(w_p)`$, describes the proportion of fish of
size $`w_p`$ in prey group $`j`$ that are hidden from predators.
$`1 - R_j(w_p)`$ is then the proportion of fish of weight $`w_p`$ and
prey group $`j`$ that are vulnerable to consumption by predators. This
model uses the method to determine the **refuge profile**, the set of
proportions that describe refuge availability across the entire size
range of model fish (Rogers et al. (2018)). The proportion of prey of
weight $`w_p`$ and group $`j`$ with access to refuge $`R_j(w_p)`$ is
given by

``` math
R_j(w_p) = min\left \{R_{max} \, , \frac{\tau\cdot\eta_k}{\sum_{i}\int_{w_{k-1}}^{w_k} N_i(w)~dw}
            ~~~~~~~ w_p \in (~w_{k-1}, w_k~] \right \}
```
{#eq-refuge_data}

The parameter $`\tau`$ is the proportion of fish with access to refuge
that are expected to utilize it, $`\eta_k`$ is the density
($`no./m^{2}`$) of refuges in size range $`(w_{k-1}, w_k]`$ and
$`\sum_{i} \int_{w_{k-1}}^{w_k} N_i(w)~dw`$ gives the total density
($`no./m^{2}`$) of fish from any group in size range $`(w_{k-1}, w_k]`$.
This represents the density of competitors for refuges in size class
$`k`$. With this method, the refuge profile is density dependent.

All fish are assumed to utilize refuge ($`\tau = 1`$). Fish smaller than
0.1 g are assumed to be larval reef fish that have not yet settled to
the reef. A maximum proportion of 98 % of fish are protected at any
given time.

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
their specific traits. The table below indicates how each species group
interacts with structures that provide predation refuge.

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

the figure below shows the density-dependent refuge profile produced by
the competitive method for the simple trait-based Bonaire model at
steady state.

![Proportion of each refuge-using group protected from predation, by
body length, at steady state, for the 10-group Karpata Reef
model.](karpata_model-description_files/figure-html/unnamed-chunk-8-1.png)

Proportion of each refuge-using group protected from predation, by body
length, at steady state, for the 10-group Karpata Reef model.

#### Consumption

Invertebrate consumption of the detrital resource and planktivory are
subject to a Holling functional response type II to represent satiation.
This relationship is defined by the maximum intake rate, which increases
allometrically with body size at rate $`n`$. The parameter $`h`$ is the
max consumption rate for an invertebrate consumer of size 1 gram. Values
for $`h`$ were chosen so that invertebrates are neither too starved nor
totally satiated.

Only a proportion $`\alpha_i`$ of consumed biomass is retained, while a
proportion $`1-\alpha_i`$ is expelled in the form of feces, which
contribute to the detritus. See the table below for the max consumption
rate for invertebrates.

No maximum consumption rate is imposed for predatory or herbivorous
groups.

|             |         h | alpha |    n |
|:------------|----------:|------:|-----:|
| pred_plank  |  13.36972 |   0.6 | 0.75 |
| parrotfish  |  63.73181 |   0.6 | 0.75 |
| farm_damsel |  14.99456 |   0.6 | 0.75 |
| herbs       |  58.95293 |   0.6 | 0.75 |
| inverts     | 235.76163 |   0.6 | 0.75 |

##### Metabolic Losses

The energy losses to metabolic needs are comprised of two components.
Standard metabolism occurs at rate $`k_{s.i}`$, scaling allometrically
with body size at rate $`p`$. The units of the coefficients $`k_{s.i}`$
are $`\text{grams}^{1-p}`$ per year. Losses due to activity and movement
occur at rate $`k_i`$ in grams per year, scaling with body size at rate
$`1`$.

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

##### Energy Invested into Reproduction

A proportion $`\psi_i(w)`$ of the energy available for growth and
reproduction is used for reproduction. This proportion changes from zero
below the weight $`w_{mat.i}`$ of maturation to one at the maximum
weight $`w_{max.i}`$, where all available energy is used for
reproduction.

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

##### Somatic growth

Energy that is left over after metabolism and reproduction is invested
in somatic growth. The growth rate of an individual of from group $`i`$
and weight $`w`$ is
``` math
  g_i(w) = E_{r.i}(w)\left(1-\psi_i(w)\right).
```
{#eq-growth} When food supply does not cover the requirements of
metabolism and activity, growth and reproduction stops. There is no
negative growth or starvation mortality.

The values for the model parameters were chosen so that the resulting
growth curves would be close to von Bertalanffy growth curves. The
parameters in the table below were taken from the literature.

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

Here the parameters $`a`$ and $`b`$ are parameters for the allometric
weight-length relationship $`w = a l^b`$ where $`w`$ is measured in
grams and $`l`$ is measured in centimetres.

### Mortality

The mortality rate $`\mu_i(w)`$ of an individual of group $`i`$ and
weight $`w`$ has three sources: predation mortality $`\mu_{p.i}(w)`$,
external mortality $`\mu_{ext.i}(w)`$, fishing mortality
$`\mu_{f.i}(w)`$. External mortality is composed of residual natural
mortality $`\mu_{nat.i}(w)`$, and senescence $`\mu_{sen.i}(w)`$. The
mortality rate for group $`i`$ is then given by:
``` math
\mu_i(w)=\mu_{p.i}(w)+\mu_{ext.i}(w)+\mu_{f.i}(w)
```
{#eq-mort}

#### Predation mortality

All consumption by fish translates into corresponding predation
mortalities on the ingested prey individuals. the rate at which all
predators from $`j`$ consume prey of size $`w_p`$ is
``` math
  \mathtt{pred\_rate}_j(w_p) = \int \phi_j(w,w_p) \gamma_j(w) N_j(w) \, dw.
```
{#eq-pred_rate}

The mortality rate due to predation is then obtained as
``` math
  \mu_{p.i}(w_p) = \sum_j \mathtt{pred\_rate}_j(w_p)\, V_{ji}(w_p)\, \theta_{ji}.
```
{#eq-mup}

#### Fishing Mortality

Like in `mizer`, fishing mortality in `mizerReef` is imposed by fishing
gears. The total per-capita fishing mortality (1/year) is obtained by
summing over the mortality from all gears,

``` math
\mu_{f.i}(w) = \sum_g F_{g,i}(w)
```

where the fishing mortality $`F_{g,i}(w)`$ imposed by gear $`g`$ on
group $`i`$ at size $`w`$ is calculated as:

``` math
F_{g,i}(w) = S_{g,i}(w) Q_{g,i} E_{g}
```

The constant $`S`$ is the selectivity by group, gear and size, $`Q`$ is
the catchability by group and gear and $`E`$ is the fishing effort by
gear.

Fishing parameters for mizerReef models can be set up with
[`setFishing()`](https://sizespectrum.org/mizer/reference/setFishing.html),
which contains the details of how to set up gears with different
selectivities and catchabilities for each group.

Fishing mortality $`\mu_{f.i}(w)`$ is calculated with the function
[`getFMort()`](https://sizespectrum.org/mizer/reference/getFMort.html).
The vignettes were run with following fishing parameters (the table
below).

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

#### External Mortality

##### Residual natural mortality

Mortality caused by illness, fishing, or predators not explicitly
included in the model is captured by $`\mu_{nat.i}(w)`$, which is
independent of group identities and abundances. It is assumed to
decrease allometrically with body size:

``` math
\mu_{nat.i}(w) = \mu_{nat} w^{1 - n}
```
{#eq-z0} where $`\mu_{nat}`$ is the residual natural mortality rate and
$`n`$ is the allometric scaling exponent. We use a residual mortality
rate of $`\mu_{nat} =`$ 0.2 per year at size 1 gram.

##### Senescence mortality

Senescence mortality $`\mu_{sen.i}(w)`$ is intended to capture mortality
due to illness or old age. It is independent of group abundances.
Senescence mortality is assumed to increase allometrically with body
size. The rate of senescence mortality is given by:

``` math
\mu_{sen.i}(w) = k_{sen} \left(\frac{log_{10}(w)}
                                    {log_{10}(w_{max.i})}\right)^{p_{sen}}
```
{#eq-extmort} with the ratio floored at zero for $`w < 1`$ g. Here
$`k_{sen}`$ (`sen_prop`) is the rate the curve approaches as $`w \to
w_{max.i}`$, and $`p_{sen}`$ (`sen_curve`) controls the steepness with
which mortality climbs as individuals approach their maximum size, and
$`w_{max.i}`$ is the maximum body size of species group $`i`$ in grams.

Senescence mortality due to illness and old age uses $`k_{sen} =`$ 0.1
and $`p_{sen} =`$ 0.3 for this model.

### Reproduction

The reproduction parameters $`\epsilon_i`$ and $`R_{max.i}`$ are not
directly observable. The values were instead chosen so as to produce
steady-state abundances of the groups that are in line with observations
and to give reasonable values for the reproduction level.

the table below gives the steady-state reproduction level which is
defined as the ratio between the actual reproduction rate $`R_i`$ and
the maximal possible reproduction rate $`R_{\max.i}`$.

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

The fish spectrum is fed by three background resources: plankton ,
algae, and detritus. Small individuals of all species will feed on
plankton while only herbivorous groups and invertebrates feed on algae
or detritus. The general mathematical form of the algae and detritus
dynamics – a linear production/consumption balance for each pool, solved
analytically each time step – was adapted from Delius et al. (2022, de
Juan et al. 2023). This is not a direct port: mizerShelf represents its
detritus/carrion pools by rescaling the total biomass of mizer’s
size-structured background resource spectrum, whereas mizerReef’s algae
and detritus are genuinely unstructured scalar pools with their own
dedicated encounter-rate matrices, added as independent components via
[`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html).
mizerReef also adds coral-reef-specific behaviour with no equivalent in
mizerShelf, most notably that algae production is treated as a fixed
rate of primary production decoupled from consumer demand rather than
tuned to match consumption (see
[`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)’s
documentation for the ecological rationale and citations).

#### Plankton

The plankton spectrum $`N_P(w)`$ tracks the abundance all planktonic
food sources. This spectrum starts at a smaller size than the fish
spectrum in order to provide food for the smallest individuals (larvae)
of the fish spectrum.

The plankton spectrum ranges from $`w_0=9\times 10^{-13}`$ to
$`w_{cutoff}=0.1`$ grams. The steady state abundance of plankton at size
1 gram is $`\kappa=11.4[g/m^{-2}]`$. The slope of the plankton spectrum
is set to $`\lambda=2.05`$.

#### Algae

The algal resource is described only by its total biomass $`B_A`$.
Feeding on algae is not is not size-based. Herbivores can feed on algae
of any size. In the steady state the total algal biomass per square
meter is $`B_A = 4.15\times 10^{-9}`$ grams.

##### Algal consumption

For each consumer species $`i`$, a parameter $`\rho_{A.i}`$ determines
the rate at which individuals of that species encounter algal biomass.
The parameters $`\rho_{i.A}`$ have units of $`g^{-n}`$ per year. They
are non-zero only for species that consume at least some algae. The
preference $`\theta_{i.A}`$ for algae is a value between 0 and 1
specifying the proportion of consumer diets that are comprised of algal
matter.

|             |         rho | interaction_algae |
|:------------|------------:|------------------:|
| parrotfish  | 40649043272 |               0.5 |
| farm_damsel |  9563740389 |               0.5 |
| herbs       | 37601006284 |               0.5 |

##### Algal production

The rate at which algal biomass is produced by the ecosystem is given by
a constant growth rate with units of grams per unit area per unit time.
At steady state algal production is 2000 grams per square meter per
year.

#### Detritus

Detritus is consumed by herbivores and benthic invertebrates. Feeding on
detritus is not size-based as fish can feed on detritus particles of any
size. The detritus resource is described only by its total biomass
$`B_D`$. In the steady state the total detrital biomass per square meter
is $`B_D = 1.143\times 10^{-10}`$ grams.

##### Detrital consumption

For each consumer species $`i`$, a parameter $`\rho_{D.i}`$ determines
the rate at which individuals of that species encounter detritus. The
parameters $`\rho_{i.D}`$ have units of $`g^{-n}`$ per year. They are
non-zero only for species that consume at least some detritus. The
preference $`\theta_{i.D}`$ for detritus is a value between 0 and 1
specifying the proportion of consumer diets comprised of detritus.

|             |          rho | interaction_detritus |
|:------------|-------------:|---------------------:|
| pred_crypt  |  34409747249 |                  0.5 |
| parrotfish  |  40649043272 |                  0.5 |
| farm_damsel |   9563740389 |                  0.5 |
| herbs       |  37601006284 |                  0.5 |
| inverts     | 300744178552 |                  1.0 |

##### Detrital production

Detritus comes from defecation and decomposing dead organisms that die
as a result of sources of external mortality. At steady state detritus
from decomposing dead organisms is 46.09 g/year and defecation produces
284.23 g/year. A proportion `ext_decomp` of the external mortality and a
proportion `sen_decomp` of senescence mortality decomposes to detritus.
The proportion of detritus that sinks from external mortality is 80 %
and the proportion that decomposes from senescence mortality is 80 %.
These values are based on estimates from Hatcher (1988).

External detritus is generated by unmodelled sources or sinks in from
the pelagic zone. At steady state external detritus production is
-305.69 g/year. This value is set so that steady state abundances match
empirical observations. Where it is negative, feces and mortality are
producing more detritus than can be consumed. Detrital biomass is
assumed to be washed away by wave action in this case.

## Tuning the Steady State

The following R script was used to tune the steady state parameters for
this model.

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

## Set model ----------------------------------------
params <- newReefParams(species_params = karpata_10plus,
                        interaction = karpata_int,
                        method = "binned",
                        method_params = tuning_profile)

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
plotSpectra(params, power = 1)

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
plotSpectra(params, total = TRUE, power = 1)
plotSpectra(params, total = TRUE, power = 2)

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
plotSpectra(params, total = TRUE, power = 2)
plotPredMort(params) + facet_wrap(~Species)
plotFeedingLevel(params)
plotDiet(params) + scale_x_log10(limits = c(1, 1e4))
plotSpectra(params, power = 1)

# Save!
caribbean_10_model <- reefSteady(params)
```

de Juan, S., G. Delius, and F. Maynou. 2023. [A model of size-spectrum
dynamics to estimate the effects of improving fisheries selectivity and
reducing discards in mediterranean mixed demersal
fisheries](https://doi.org/10.1016/j.fishres.2023.106764). Fisheries
Research 266:106764.

Delius, G., S. de Juan, and F. Maynou. 2022. [mizerShelf: Mizer models
with carrion and detritus components suitable for continental shelf
ecosystems](https://sizespectrum.org/mizerShelf/).

Delius, G., F. Scott, J. Blanchard, and K. Andersen. 2023. [Mizer:
Dynamic multi-species size spectrum
modelling](https://CRAN.R-project.org/package=mizer).

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
