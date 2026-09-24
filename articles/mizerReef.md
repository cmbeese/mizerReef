# Getting started with MizerReef

## Overview

The mizerReef package enables multi-species dynamic size-spectrum
modelling in R, with an explicit, mechanistic representation of habitat
structural complexity. In this vignette, we walk through the basic steps
needed to build, tune and explore your first MizerReef model, including:

1.  [Installing MizerReef](#installing-mizerreef)
2.  [Setting species parameters](#setting-species-parameters)
3.  [Setting the refuge profile](#setting-the-refuge-profile)
4.  [Creating your first model](#creating-your-first-model)
5.  [Tuning the steady state](#tuning-the-steady-state)
6.  [Exploring results](#exploring-results)

Once you have a tuned model,
[`vignette("running-simulations")`](https://cmbeese.github.io/mizerReef/articles/running-simulations.md)
picks up from there: changing the refuge profile and projecting a model
forward in time (fishing pressure, habitat degradation trajectories, and
so on).

> **Model context:** Habitat structure mediates system dynamics by
> providing predation refuge. The most effective refuges fit prey while
> excluding predators, so refuge use is primarily governed by body size.
> In systems with high structural complexity (for example, coral reefs),
> these patterns are reflected in the size structure of fish
> assemblages.
>
> Because refuge protection is size-dependent, size-spectrum models
> provide a natural framework for exploring how benthic structure
> influences community dynamics. MizerReef modifies predator–prey
> encounter rates to represent the effects of habitat structure,
> allowing users to explicitly account for changes in refuge
> availability caused by habitat degradation or modification.

For a detailed description of the model formulation and supporting
references, see Chapter 3 of [Modelling Coral Reef Futures: Exploring
the role of structural complexity in sustaining ecosystem services (PhD
Thesis,
VUW)](https://openaccess.wgtn.ac.nz/articles/thesis/Modelling_Coral_Reef_Futures_Exploring_the_role_of_structural_complexity_in_sustaining_ecosystem_services/26421523?file=48064144).

MizerReef builds on the `mizer` package

`MizerReef` is an extension of the `mizer` package and uses many of the
same functions and parameters. If you are new to `mizer`, want a
refresher on the general workflow, or would like more background on
size-spectrum modelling, visit the [mizer
website](https://sizespectrum.org/mizer/). You may also find the [mizer
course](https://mizer.course.nov22.sizespectrum.org/) helpful.

## Installing MizerReef

[Skip to Setting up species parameters](#setting-up-species-parameters)

MizerReef is currently only available from GitHub. To install the latest
version, use the `devtools` package:

``` r

install.packages("devtools")
devtools::install_github("cmbeese/mizerReef")
```

MizerReef depends on the `mizer` and `mizerExperimental` packages. If
not already installed, R will prompt you to install them automatically.
Without them, you will not be able to use all of the features in
`mizerReef`.

After installation, load Mizer, mizerExperimental, and MizerReef in each
new R session:

``` r

library(mizer)
library(mizerExperimental)
library(mizerReef)
```

> 💡 **Tip:** Be sure to **load mizerReef last**. Some of the functions
> in mizerReef override functions in mizer, so loading in this order
> ensures that correct versions are used.

Use a recent version of R and RStudio for best results. For
troubleshooting or more details, see the [MizerReef
documentation](https://cmbeese.github.io/mizerReef/) or the [GitHub
repository](https://github.com/cmbeese/mizerReef).

## Setting species parameters

[Skip to mizerReef
parameters](#additional-species-parameters-needed-for-mizerreef)

The species parameter table is the foundation of any mizer or mizerReef
model. It describes the biological and ecological traits of each species
in your system.

> 💡 **Tip:** Many users prefer to create this table in a spreadsheet
> program (like Excel or Google Sheets) and then import it into R as a
> data frame.

### Species parameters required by base mizer

For a multi-species mizer model, only the following columns are strictly
required in the species parameter data frame:

- `species`: Name of the species or group
- `w_max` or `l_max`: Maximum observed weight or length (if providing
  length, you must also provide length-weight conversion parameters `a`
  and `b`)

By default, mizerReef creates multispecies mizer models.

mizerReef is not compatible with trait-based or community models at this
time. See [multispecies mizer
models](https://sizespectrum.org/mizer/articles/multispecies_model.html#overview)
to learn more about the differences between the three model types in
mizer.

You can tune your model to a specific system using abundance data if you
also provide:

- `biomass_observed`: Observed abundance for each species.
- `biomass_cutoff`: Minimum weight of organisms caught by survey methods
  (helps with tuning)

The choice in units for your data is arbitrary as long as you are
consistent.

Abundance data can be given as numbers per area, numbers per volume or
total numbers for the entire study area. See [Units in
Mizer](https://sizespectrum.org/mizer/reference/setParams.html#units-in-mizer)
for more information.

Since MizerReef’s vulnerability and refuge dynamics depend on size, it
is good to include parameters related to growth and size rather than
relying on defaults, including:

- `w_mat`: Maturity weight (important for life history & growth)
- `beta` & `sigma`: Lognormal predation kernel parameters (set for each
  species, can use other kernels)
- Length-to-weight conversion parameters `a` and `b`

You should also provide the interaction matrix, which specifies
predator-prey relationships between species.

Each value in the [interaction
matrix](https://sizespectrum.org/mizer/reference/setParams.html#setting-interaction-matrix),
ranging from 0 to 1, represents the strength of interaction between a
predator (row) and its prey (column). These can represent spatial
overlap, diet preferences, or other ecological factors influencing
predation rates.

It is important that the order of rows and columns in the interaction
matrix matches the order of species in the species parameter data frame
within the params object. To view the [example interaction
matrix](https://cmbeese.github.io/mizerReef/reference/caribbean_3_interaction.html)
included with MizerReef, run:

``` r

data("caribbean_3_interaction")
caribbean_3_interaction
```

See [Setting Parameters in
Mizer](https://sizespectrum.org/mizer/reference/setParams.html) for
details on additional optional columns and their defaults.

### Additional species parameters needed for mizerReef

[Skip to Example species parameters](#example-species-parameters)

mizerReef extends mizer by modeling reef-specific dynamics. To utilise
the predation vulnerability and unstructured resource dynamics in
MizerReef, add these columns to your species parameter data frame:

| Column | Type | Description |
|:---|:---|:---|
| refuge_user | logical | TRUE if the group uses predation refuge (i.e., individuals hide inside habitat structure to avoid predators; typical for small-bodied or cryptic reef fish). FALSE for species that do not use refuge. FALSE by default |
| blocked_pred | logical | TRUE if the group is blocked from accessing prey in refuge (i.e., predators that cannot reach prey hiding in refuge). FALSE for species with behavioural or morphological adaptations (e.g., eels) that allow them to access prey in refuge. FALSE by default. |
| satiation | logical | TRUE if group is subject to satiation on unstructured resources (algae or detritus). By default TRUE only for detritivores – species that eat detritus but do not also graze algae or eat other species – and FALSE for every other group, including herbivores. |
| interaction_algae | numeric | Proportion of diet from algae (0–1). 0 by default. |
| interaction_detritus | numeric | Proportion of diet from detritus (0–1). 0 by default |

Additional columns required for mizerReef species parameters. {.table}

Every one of these is a **mizerReef default**, not a mizer one: leave a
column out (or leave a cell blank/`NA`) and mizerReef fills it in with
the value above, switching that feature off for the species and
reporting that it did so. Override any of them by including the column
in your own species table, or afterwards with
`species_params(params)$refuge_user <- ...` (and similarly for the
others) or the relevant setter
([`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)
for the first three,
[`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)/[`setDetritusParams()`](https://cmbeese.github.io/mizerReef/reference/setDetritusParams.md)
for the interaction columns). Those two setters also fill in several
*model-level* defaults – not tied to one species, e.g. how much of a
refuge is ever protected, or algae’s production rate – covered next in
[Setting the refuge profile](#setting-the-refuge-profile) and in
[`vignette("tuning-diet-composition")`](https://cmbeese.github.io/mizerReef/articles/tuning-diet-composition.md).

> 💡 **Tip:** When importing your species parameter table from a CSV,
> make sure that missing values are represented as NA or blank cells,
> not zeros.
>
> You can load your species parameter data frame into R using standard
> functions like [`read.csv()`](https://rdrr.io/r/utils/read.table.html)
> or `readxl::read_excel()`, depending on your file format. Ensure that
> the column names in your data frame match those expected by MizerReef.
> Use `na.strings = c(““,” “)` in
> [`read.csv()`](https://rdrr.io/r/utils/read.table.html) to ensure
> blanks are read as `NA` rather than `0`.

### Example mizerReef species parameters

[Skip to Setting the refuge profile](#setting-the-refuge-profile)

MizerReef includes example model parameters with the package. See the
[Model
Description](https://cmbeese.github.io/mizerReef/articles/karpata_model-description.md)
for more details on these example data.

The included species parameter data frame (`caribbean_10_species`) is
based on fish assemblage data from a Caribbean reef site with relatively
low fishing pressure. The full species parameter data frame for the
Caribbean 10 example is shown below:

| species | l_max | w_mat | age_mat | beta | sigma | biomass_cutoff | biomass_observed | a | b | interaction_detritus | interaction_algae | refuge_user | blocked_pred | satiation |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|:---|:---|:---|
| pred_eng | 45 | 50.0 | 2.0 | 60 | 2 | 12.62969 | 5.51 | 0.01100 | 3.06 | 0.0 | 0.0 | TRUE | TRUE | FALSE |
| pred_grab | 60 | 140.0 | 5.4 | 40 | 2 | 17.80530 | 27.60 | 0.01740 | 3.01 | 0.0 | 0.0 | TRUE | TRUE | FALSE |
| eels | 100 | 300.0 | 3.0 | 30 | 2 | NA | 10.00 | 0.00098 | 3.24 | 0.0 | 0.0 | TRUE | FALSE | FALSE |
| pred_crypt | 8 | 0.8 | NA | 10 | 1 | 1.00000 | NA | 0.01122 | 3.04 | 0.5 | 0.0 | TRUE | FALSE | FALSE |
| pred_inv | 45 | 50.0 | 0.5 | 50 | 2 | 26.58539 | 14.13 | 0.01200 | 3.10 | 0.0 | 0.0 | TRUE | FALSE | FALSE |
| pred_plank | 20 | 2.5 | 1.0 | 1000 | 3 | 13.49043 | 0.55 | 0.01259 | 3.03 | 0.0 | 0.0 | TRUE | FALSE | TRUE |
| parrotfish | 64 | 63.0 | 1.6 | 30 | 1 | 15.48385 | 30.56 | 0.01380 | 3.05 | 0.5 | 0.5 | TRUE | FALSE | TRUE |
| farm_damsel | 13 | 1.0 | 1.0 | 30 | 1 | NA | 0.40 | 0.02042 | 2.97 | 0.5 | 0.5 | TRUE | FALSE | TRUE |
| herbs | 39 | 105.0 | 2.0 | 30 | 1 | 75.75344 | 1.50 | 0.02570 | 2.95 | 0.5 | 0.5 | TRUE | FALSE | TRUE |
| inverts | 30 | 0.1 | NA | 30 | 2 | 3.12500 | NA | 0.02500 | 3.00 | 1.0 | 0.0 | FALSE | FALSE | TRUE |

> 💡 **Tip:** The Karpata species parameters also include optional
> columns like `w_mat`, `age_mat`, `k_vb`, and `ks`. Include as many
> parameters as you have data for to assist in the calibration process.

## Setting the refuge profile

[Skip to Creating your first model](#creating-your-first-model)

The refuge profile defines how predation refuge availability varies with
prey size (see Figure 1). This is a key feature of MizerReef that allows
you to represent the effects of habitat structure on predator-prey
interactions.

![Schematic plot: x axis is log body size, y axis is proportion
protected. Red fish icons show fish sizes, with many small fish and one
large fish. Diagonal line shows decreasing protection with size.
Transparent grey bars across size bins indicate proportion protected,
representing protection across the predicted Sheldon
spectrum.](figures/refuge_profile.png)

Figure 1. Conceptual schematic of the refuge profile. The red fish
represents the modelled fish spectrum. Individuals protected from
predators by refuge are covered by the grey box. The remaining
individuals are vulnerable to predation.

`MizerReef` currently provides three methods to define how refuge
availability varies with prey size:

**Sigmoidal**: Good for data-poor reefs or when you want a simple,
smooth profile.

- a smooth declining function controlled by a threshold length
  (L_refuge) and a maximum proportion protected

  - `method_params` should be a **list** or **data frame** with:
    - `L_refuge`: numeric, threshold length (cm) at which refuge
      protection starts to decline.
    - `prop_protect`: numeric, maximum proportion of individuals
      protected by refuge (0–1).

  **Example:**

  ``` r

  method_params = list(L_refuge = 10, prop_protect = 0.8)
  ```

**Binned**: Good for theoretical experiments or when you have coarse bin
information.

- user-specified length bins with a constant protection proportion
  inside each bin

  - `method_params` should be a **data frame** or **matrix** with two
    columns:
    - `length_bin`: numeric, the upper length (cm) of each bin.
    - `protection`: numeric, the proportion of individuals protected by
      refuge (0–1) within each length bin.

  **Example:**

  ``` r

  method_params = data.frame(
    length_bin = c(5, 10, 20, 40),
    protection = c(1, 0.5, 0, 0.2)
  )
  ```

**Competitive**: divides refuges among similarly sized competitors. Use
this when you have empirical refuge density data.

- uses refuge density (no./m^2) for each length bin, protection depends
  on fish density within each bin (density-dependent)

  - `method_params` should be a **data frame** or **matrix** with two
    columns:
    - `length_bin`: numeric, the upper length (cm) of each bin.
    - `refuge_density`: numeric, the density of refuges (no./m^2)
      available for each length bin.

  **Example:**

  ``` r

  method_params = data.frame(
    length_bin = c(5, 10, 20, 40),
    refuge_density = c(2, 1, 0.5, 3)
  )
  ```

A few more refuge settings apply across all three methods, and mizerReef
fills them in with a default if you don’t set them via
[`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)/[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md):

| Argument | Default | Description |
|:---|:---|:---|
| max_protect | 0.98 | Cap on the proportion of any size class that refuge can ever protect, so some food always remains accessible to predators. |
| tau | 1 | Proportion of refuge-eligible individuals that actually use a refuge when one is available. |
| w_settle | 0.1 g | Weight below which fish are assumed unsettled larvae, exempt from the refuge calculation. |
| a_bar / b_bar | 0.025 / 3 | Length-weight parameters for a ‘dummy fish’, used to convert refuge bin boundaries between length and weight, and to fill in a real species’ own a/b if it’s missing one. |
| use_dummy_fish_bins | TRUE | Whether refuge bins are set by weight, using a_bar/b_bar (TRUE), or by length, using each species’ own a/b (FALSE). |

mizerReef’s default refuge settings, shared across all three methods.
{.table}

`a_bar`/`b_bar` and `w_settle` were chosen for small-bodied coral-reef
fish, not for fish in general – a pelagic species’ length-weight
relationship, for instance, can be very different from `a_bar = 0.025`,
`b_bar = 3`. If your system isn’t a coral reef, treat every default on
this page as a starting point to check against your own species’ traits,
not a universal constant.

### Refuge profiles and body shape

This plot shows example refuge profiles created using each method and
how they differ based on species body shape characteristics. The same
set of species groups can receive different protection based on their
body shape, the chosen method, and method-specific parameters.

![Four-panel plot: each panel shows protection proportion across body
size for a different species with distinct body shape, comparing three
refuge profile methods. Panels illustrate how protection varies for
deep, compressed, elongate, and fusiform
species.](figures/body-shape-example.png)

Figure 2. Example refuge profiles for three methods (sigmoidal, binned,
competitive) applied to species with different body shapes (deep,
compressed, elongate, fusiform).

The package includes several example refuge profiles for tuning and
demonstration.

mizerReef’s example models use competitive refuge profiles based on
field data. Use the code below to view the built-in Karpata Reef refuge
profile:

``` r

data(karpata_refuge)
karpata_refuge
```

    ##    start_L end_L refuge_density
    ## 1        0     5     7.53333333
    ## 2        5    10     1.40000000
    ## 3       10    15     0.70833333
    ## 4       15    20     0.28333333
    ## 5       20    25     0.10000000
    ## 6       25    30     0.05000000
    ## 7       30    35     0.04166667
    ## 8       35    40     0.03333333
    ## 9       40    45     0.03333333
    ## 10      45    50     0.04166667

See [example
models](https://cmbeese.github.io/mizerReef/reference/index.html#example-models)
for more details on built-in refuge profiles.

## Creating your first model

[Skip to Tuning the steady state](#tuning-the-steady-state)

Once you have your species parameters (with reef-specific columns) and
interaction matrix ready, you can create a `MizerParams` object using
the
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
function. This function extends `mizer`’s
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.html)
by adding reef-specific arguments, checking user-supplied parameters,
and setting sensible defaults for any missing reef-specific values.

> 💡 **Tip:** When creating a new `mizerReef` model with
> [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md),
> you can’t use the competitive refuge method when calibrating biomasses
> because it is density-dependent.

Use a tuning profile instead.

The best practice if you have refuge density data is to first create a
model using the binned method that approximates your refuge profile to
reach an initial steady state and calibrate biomasses, then switch to
the competitive method. The package includes an example refuge profile
for tuning (`tuning_profile`) that can be used for this purpose.

``` r

caribbean_10_model <- newReefParams(species_params = caribbean_10_species,
                                    interaction = caribbean_10_interaction,
                                    method = "binned",
                                    method_params = tuning_profile)
```

Besides the species-level and refuge-level defaults already covered
above,
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
(all three of these are set directly in its own function signature/body,
not inherited from a setter) makes a few default choices of its own
about the underlying mizer model:

| Argument | mizerReef default | mizer’s own default | Description |
|:---|:---|:---|:---|
| w_pp_cutoff | 1 g | 10 g | Maximum size of the plankton resource spectrum. |
| n (and p) | 0.75 | 2/3 (≈0.667); p defaults independently to 0.7 | Allometric growth exponent; mizerReef also uses it for the metabolic exponent p, rather than defaulting p independently. |
| crit_feed | 0.6 | 0.6, for the f0 species parameter | Target feeding level; fills in f0 for any species that doesn’t have its own, which in turn calibrates that species’ algae/detritus encounter rate (rho) at maximum size. |

mizerReef’s own defaults for underlying mizer construction arguments.
{.table}

These are ordinary arguments, so pass your own value the same way as any
other (e.g. `newReefParams(..., w_pp_cutoff = 10)`). They don’t all have
the same kind of justification, though:

- **`w_pp_cutoff` is a deliberate, architectural choice**, not an
  oversight. mizerReef’s narrower plankton spectrum represents plankton
  *only* – unlike a typical mizer model, where the single background
  resource often also stands in for small invertebrates, a reef model
  gives invertebrates their own explicit species/spectrum instead, so
  the plankton resource doesn’t need to reach as high in size.
- **`crit_feed` matches mizer’s own `f0` default** (it used to be `0.7`,
  with no documented reason for the divergence – fixed in 2.0.3).
- **`n`/`p` has no documented rationale for differing from mizer’s own
  default.** It traces to the trait-based model in Rogers (2018) (see
  [`vignette("caribbean_3_model-description")`](https://cmbeese.github.io/mizerReef/articles/caribbean_3_model-description.md))
  rather than a reef-specific finding, so it’s a starting point worth
  sensitivity-testing for your own system, not a value to trust just
  because it’s what mizerReef ships.

For every other construction argument (`species_params`, `interaction`,
`kappa`, `resource_rate`, and so on),
[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
simply passes your value straight through to
[`newMultispeciesParams()`](https://sizespectrum.org/mizer/reference/newMultispeciesParams.html),
so mizer’s own defaults and the messages it reports when using them
apply unchanged – see mizer’s [Setting
Parameters](https://sizespectrum.org/mizer/reference/setParams.html)
reference page for those.

After creating your initial `params` object, you will typically run
through a tuning sequence to calibrate biomasses and adjust
reproduction, growth, and unstructured resource parameters to match
observed data.

## Tuning the steady state

[Skip to Exploring results](#exploring-results)

Reaching a steady state that matches observed biomasses and growth rates
is nontrivial and often unique for each system. The procedure developed
here was suitable for the Karpata reef data but may differ depending on
your calibration data. In brief, the tuning procedure is as follows:

1.  Start with plausible species parameters. Create an initial `params`
    object with
    [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
    using a binned or sigmoidal refuge profile that mimics your data.
2.  Reduce the density-dependence of reproduction by reducing the
    reproduction level. Run to a steady using reefSteady(). Check the
    resource abundance and scale if needed.
3.  Iterate through
    [`calibrateReefBiomass()`](https://cmbeese.github.io/mizerReef/reference/calibrateReefBiomass.md),
    [`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.html),
    [`matchReefGrowth()`](https://cmbeese.github.io/mizerReef/reference/matchReefGrowth.md)
    and
    [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
    to reach a satisfactory steady state.  
4.  Change to your desired refuge method (for example, competitive)
    using
    [`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md),
    then re-tune the steady state by iterating/repeating step 3.
5.  Tune the reproduction parameters according to the mizer blog recipe
    to reach the final steady state.

> 💡 **Tip:** One of the most useful summary plots for steady state
> calibration is
> [`plotBiomassObservedVsModel()`](https://sizespectrum.org/mizer/reference/plotBiomassObservedVsModel.html),
> which shows the total biomass for each species in your model
> vs. observed values.

> 💡 **Tip:** There are many reasonable ways to reach a suitable steady
> state. This recipe only represents one approach. You may need to
> adjust the steps or order depending on your system and data.
>
> The tuning approach used for MizerReef is adapted from the [5-step
> recipe](https://blog.mizer.sizespectrum.org/posts/2021-08-20-a-5-step-recipe-for-tuning-the-model-steady-state/)
> described in the mizer blog.

For more information on tuning MizerReef models, see the [MizerReef
Steady State
recipe](https://cmbeese.github.io/mizerReef/articles/steady-state-recipe.md).

## Exploring results

[Skip to Next steps](#next-steps)

After reaching a steady state, you should explore the results to ensure
they make ecological sense and match expectations for your system.
MizerReef provides several plotting functions to help you visualize and
interpret your model outputs:

- **Refuge profile at steady state:**  
  Use \[plotRefugeProfile()\] to see the proportion of individuals
  protected by refuge across sizes and species.

  ``` r

  plotRefugeProfile(params)
  ```

- **Diet composition:**  
  Use \[plotDiet()\] to display the proportion of each prey type in the
  diet of each predator.

  ``` r

  plotDiet(params)
  ```

- **Total biomass for each species:** Use \[plotBiomass()\] on a
  projected simulation to view biomasses over time.

  ``` r

  sim <- project(params)
  plotBiomass(sim)
  ```

> 💡 **Tip:** `plotBiomass()` is mizer’s own function.
>
> mizerReef does not override it — instead it registers a
> [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.html)
> method for `mizerReef` models so that algae and detritus biomass are
> automatically included alongside species biomass, regardless of the
> order in which `mizer` and `mizerReef` are loaded.

- **Productivity by species:**  
  Use \[plotProductivity()\] to view total productivity for each
  species.

  ``` r

  plotProductivity(params)
  ```

For a full list of available summary and diagnostic plots, see
[MizerReef summary
plots](https://cmbeese.github.io/mizerReef/reference/index.html#summary-plots)
and [mizer’s plotting results reference
page](https://sizespectrum.org/mizer/reference/index.html#plotting-results).

## Next steps

Once your model is calibrated and validated as above, see [Running
MizerReef
Simulations](https://cmbeese.github.io/mizerReef/articles/running-simulations.md)
for how to change the refuge profile and project a tuned model forward
in time – fishing pressure, habitat degradation trajectories, and more –
including that vignette’s own [Links and further
reading](https://cmbeese.github.io/mizerReef/articles/running-simulations.html#links-and-further-reading)
section.

## Session info

    ## R version 4.6.1 (2026-06-24)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.5 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] knitr_1.52              mizerReef_2.1.0         mizerExperimental_3.3.1
    ## [4] mizer_3.4.0.9000       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] plotly_4.12.1       sass_0.4.10         generics_0.1.4     
    ##  [4] tidyr_1.3.2         stringi_1.8.9       digest_0.6.39      
    ##  [7] magrittr_2.0.5      timechange_0.4.0    evaluate_1.0.5     
    ## [10] grid_4.6.1          RColorBrewer_1.1-3  fastmap_1.2.0      
    ## [13] plyr_1.8.9          jsonlite_2.0.0      httr_1.4.9         
    ## [16] purrr_1.2.2         viridisLite_0.4.3   scales_1.4.0       
    ## [19] textshaping_1.0.5   jquerylib_0.1.4     cli_3.6.6          
    ## [22] rlang_1.3.0         cachem_1.1.0        yaml_2.3.12        
    ## [25] otel_0.2.0          tools_4.6.1         reshape2_1.4.5     
    ## [28] dplyr_1.2.1         ggplot2_4.0.3       assertthat_0.2.1   
    ## [31] vctrs_0.7.3         R6_2.6.1            lubridate_1.9.5    
    ## [34] lifecycle_1.0.5     stringr_1.6.0       fs_2.1.0           
    ## [37] htmlwidgets_1.6.4   ragg_1.5.2          pkgconfig_2.0.3    
    ## [40] desc_1.4.3          pkgdown_2.2.1       pillar_1.11.1      
    ## [43] bslib_0.12.0        gtable_0.3.6        glue_1.8.1         
    ## [46] data.table_1.18.6.1 Rcpp_1.1.2          systemfonts_1.3.2  
    ## [49] xfun_0.61           tibble_3.3.1        tidyselect_1.2.1   
    ## [52] farver_2.1.2        htmltools_0.5.9     rmarkdown_2.32     
    ## [55] compiler_4.6.1      S7_0.2.2
