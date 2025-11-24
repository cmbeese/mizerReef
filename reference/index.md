# Package index

## Setting up a mizerReef model

These functions are for setting up a new mizerReef model and finding a
steady state for the dynamical system.

- [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
  : Set up parameters for a mizerReef model
- [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  : Project a mizerReef model to steady state

## Tuning a mizerReef model

These functions are for calibrating mizerReef models to match empirical
observations.

- [`calibrateReefBiomass()`](https://cmbeese.github.io/mizerReef/reference/calibrateReefBiomass.md)
  : Calibrate the scale of a mizerReef model to match total observed
  biomass
- [`calibrateReefNumber()`](https://cmbeese.github.io/mizerReef/reference/calibrateReefNumber.md)
  : Calibrate the model scale to match total observed number
- [`matchReefGrowth()`](https://cmbeese.github.io/mizerReef/reference/matchReefGrowth.md)
  : Match observed growth rates
- [`scaleReefAbundance()`](https://cmbeese.github.io/mizerReef/reference/scaleReefAbundance.md)
  : Scale reef abundances
- [`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md)
  : Scale model parameters
- [`step_tune`](https://cmbeese.github.io/mizerReef/reference/step_tune.md)
  : Stepped refuge profile for tuning steady states
- [`tuning_profile`](https://cmbeese.github.io/mizerReef/reference/tuning_profile.md)
  : Constant refuge profile for tuning steady states

### Tuning Profiles

These refuge profiles are useful for tuning mizerReef models to steady
state.

- [`tuning_profile`](https://cmbeese.github.io/mizerReef/reference/tuning_profile.md)
  : Constant refuge profile for tuning steady states
- [`step_tune`](https://cmbeese.github.io/mizerReef/reference/step_tune.md)
  : Stepped refuge profile for tuning steady states

## Predation Refuge

These functions set up the mediation of predation with refuge.

### Setting the refuge parameters

These functions allow users to set or change the refuge profile.

- [`getRefuge()`](https://cmbeese.github.io/mizerReef/reference/getRefuge.md)
  : Finds the refuge length bins by functional group and stores them
  params
- [`newRefuge()`](https://cmbeese.github.io/mizerReef/reference/newRefuge.md)
  : Change the refuge parameters for a model
- [`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)
  : Set the refuge profile parameters

### New Rate Functions

These functions supplement or replace mizer’s default rate functions so
that refuge impacts simulations.

- [`getVulnerable()`](https://cmbeese.github.io/mizerReef/reference/getVulnerable.md)
  : Get vulnerability level at in time range t
- [`reefEncounter()`](https://cmbeese.github.io/mizerReef/reference/reefEncounter.md)
  : Get encounter rate needed to project a mizerReef model
- [`reefPredMort()`](https://cmbeese.github.io/mizerReef/reference/reefPredMort.md)
  : Get total predation mortality rate needed to project mizer reef
  model
- [`reefRates()`](https://cmbeese.github.io/mizerReef/reference/reefRates.md)
  : Get all rates needed to project a mizerReef model
- [`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
  : Find the proportion of fish vulnerable to being encountered by
  predators at each time step

### Degradation

These functions prepare a mizer model for projections with degradation.

- [`algae_scale`](https://cmbeese.github.io/mizerReef/reference/algae_scale.md)
  : Algae trajectory refuge density scaling parameters
- [`constant_scale`](https://cmbeese.github.io/mizerReef/reference/constant_scale.md)
  : Trajectory with no refuge density scaling for testing
- [`getDegrade()`](https://cmbeese.github.io/mizerReef/reference/getDegrade.md)
  : Get vulnerability level at in time range t
- [`recovery_scale`](https://cmbeese.github.io/mizerReef/reference/recovery_scale.md)
  : Recovery trajectory refuge density scaling parameters
- [`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
  : Scales the refuge density by a given value at set times
- [`rubble_scale`](https://cmbeese.github.io/mizerReef/reference/rubble_scale.md)
  : Rubble trajectory refuge density scaling parameters
- [`setDegradation()`](https://cmbeese.github.io/mizerReef/reference/setDegradation.md)
  : Prepare a steady state model for projections with degradation

### Plotting the refuge profile

These functions allow users visualize the refuge profile for different
functional groups in terms of their body length.

- [`plotRefuge()`](https://cmbeese.github.io/mizerReef/reference/plotRefuge.md)
  [`plotlyRefuge()`](https://cmbeese.github.io/mizerReef/reference/plotRefuge.md)
  : Plot the refuge profile, species by length
- [`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
  [`plotlyVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
  : Plot the vulnerability to predation of species by weight

## Unstructured Resources

These functions set up the production and consumption dynamics for the
algae and detritus resources and provide some plotting abilities to
assess these dynamics.

### Algae

- [`algae_biomass()`](https://cmbeese.github.io/mizerReef/reference/algae_biomass.md)
  : algae Biomass
- [`algae_consumption()`](https://cmbeese.github.io/mizerReef/reference/algae_consumption.md)
  : Mass-specific algae consumption rate
- [`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md)
  : Algae dynamics
- [`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md)
  : Algae dynamics with carrying capacity
- [`getAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeConsumption.md)
  : Get algae consumption rates
- [`getAlgaeProduction()`](https://cmbeese.github.io/mizerReef/reference/getAlgaeProduction.md)
  : Algae production rate
- [`plotAlgaeConsumption()`](https://cmbeese.github.io/mizerReef/reference/plotAlgaeConsumption.md)
  : Plot algae consumption rates
- [`rescale_algae()`](https://cmbeese.github.io/mizerReef/reference/rescale_algae.md)
  : Rescale algae biomass without changing anything else

### Detritus

- [`detritus_biomass()`](https://cmbeese.github.io/mizerReef/reference/detritus_biomass.md)
  : Detritus Biomass
- [`detritus_consumption()`](https://cmbeese.github.io/mizerReef/reference/detritus_consumption.md)
  : Mass-specific detritus consumption rate
- [`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md)
  : Detritus dynamics
- [`detritus_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics_cc.md)
  : Detritus dynamics with carrying capacity
- [`detritus_lifetime()`](https://cmbeese.github.io/mizerReef/reference/detritus_lifetime.md)
  [`` `detritus_lifetime<-`() ``](https://cmbeese.github.io/mizerReef/reference/detritus_lifetime.md)
  : Expected detritus lifetime
- [`getDetritusConsumption()`](https://cmbeese.github.io/mizerReef/reference/getDetritusConsumption.md)
  : Get detritus consumption rates
- [`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md)
  : Detritus production rate
- [`plotDetritusConsumption()`](https://cmbeese.github.io/mizerReef/reference/plotDetritusConsumption.md)
  : Plot detritus consumption rates by species
- [`plotDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/plotDetritusProduction.md)
  : Plot detritus production rates from each source
- [`rescale_detritus()`](https://cmbeese.github.io/mizerReef/reference/rescale_detritus.md)
  : Rescale detritus biomass without changing anything else

### Other Components

- [`constant_dynamics()`](https://cmbeese.github.io/mizerReef/reference/constant_dynamics.md)
  : Hold resource dynamics constant
- [`encounter_contribution()`](https://cmbeese.github.io/mizerReef/reference/encounter_contribution.md)
  : Contribution of unstructured components to the encounter rate
- [`rescaleComponents()`](https://cmbeese.github.io/mizerReef/reference/rescaleComponents.md)
  : Rescale algae and detritus biomass without changing anything else
- [`scaleReefBackground()`](https://cmbeese.github.io/mizerReef/reference/scaleReefBackground.md)
  : Scale background down by a factor
- [`setURParams()`](https://cmbeese.github.io/mizerReef/reference/setURParams.md)
  : Checks unstructured resource parameters and interaction matrix
- [`setURcapacity()`](https://cmbeese.github.io/mizerReef/reference/setURcapacity.md)
  : Switch to unstructured resource dynamics with carrying capacities
- [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)
  : Tune unstructured resources (algae and detritus) to steady state
- [`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
  : Tune unstructured resources with carrying capacities (algae and
  detritus) to steady state

## External Mortality

These functions add size-dependent senescence mortality to mizerReef
models.

- [`getSenMort()`](https://cmbeese.github.io/mizerReef/reference/getSenMort.md)
  : Get the size specific senescence mortality rate
- [`reefFeedingLevel()`](https://cmbeese.github.io/mizerReef/reference/reefFeedingLevel.md)
  : Reef feeding level
- [`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md)
  : Total mortality rate in the reef ecosystem model
- [`reefSenMort()`](https://cmbeese.github.io/mizerReef/reference/reefSenMort.md)
  : Expanding external mortality rate to include senescence
- [`setExtMortParams()`](https://cmbeese.github.io/mizerReef/reference/setExtMortParams.md)
  : Set the parameters for external mortality

## Summary Functions

- [`getEGrowthTime()`](https://cmbeese.github.io/mizerReef/reference/getEGrowthTime.md)
  : Get energy rate available for growth through time
- [`getProductivity()`](https://cmbeese.github.io/mizerReef/reference/getProductivity.md)
  : Calculate fisheries productivity for each species group

## Summary Plots

These functions calculate and plot summary statistics and allow for the
comparison of results between different models.

- [`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md)
  [`plotly2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md)
  : Plot the fisheries productivity of two models or two different size
  ranges in the same plot
- [`plot2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md)
  [`plotly2TotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plot2TotalBiomass.md)
  : Plot the total biomass of two models or of two different size ranges
  in the same plot
- [`plotBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotBiomass.md)
  [`plotlyBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotBiomass.md)
  : Plot the biomass of Species Groups and unstructured components
  through time
- [`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)
  [`plotlyProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)
  : Plot the total productivity for each species Group
- [`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md)
  [`plotlyProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md)
  : Plot the relative difference between the potential fisheries
  productivity rates of two models or two different size ranges in the
  same plot
- [`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md)
  : Plot the relative contribution of each species group to total
  abundance, total biomass, and total productivity
- [`plotSpectra2()`](https://cmbeese.github.io/mizerReef/reference/plotSpectra2.md)
  : Show two size spectra in the same plot
- [`plotSpectraRelative()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraRelative.md)
  [`plotlySpectraRelative()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraRelative.md)
  : Plot the relative difference or percent change between two spectra
- [`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md)
  : Plot the total abundance for each species group at steady state
- [`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)
  [`plotlyTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)
  : Plot the total fishable biomass for each Species Group at steady
  state
- [`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md)
  [`plotlyTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md)
  : Plot the relative difference in between the total fishable biomasses
  of each each Species Group at steady state

## Helper Functions

- [`different()`](https://cmbeese.github.io/mizerReef/reference/different.md)
  : Check whether two objects are different
- [`removeSpecies()`](https://cmbeese.github.io/mizerReef/reference/removeSpecies.md)
  : Remove some species from the model

## Example Models

These files hold example parameters and mizer params object to explore
models set up to emulate the coral reefs of Bonaire.

### Simple trait based model

- [`bonaire_int`](https://cmbeese.github.io/mizerReef/reference/bonaire_int.md)
  : interaction matrix for for a simple, generic Caribbean reef
- [`bonaire_refuge`](https://cmbeese.github.io/mizerReef/reference/bonaire_refuge.md)
  : Competitive method refuge parameters for a simple, generic Caribbean
  reef
- [`bonaire_species`](https://cmbeese.github.io/mizerReef/reference/bonaire_species.md)
  : species_params dataframe for a simple, generic Caribbean reef
- [`bonaire_model`](https://cmbeese.github.io/mizerReef/reference/bonaire_model.md)
  : MizerParams object for a simple, generic Caribbean reef

### Multi-species model

- [`karpata_int`](https://cmbeese.github.io/mizerReef/reference/karpata_int.md)
  : Interaction matrix for a generic Caribbean reef
- [`karpata_refuge`](https://cmbeese.github.io/mizerReef/reference/karpata_refuge.md)
  : Competitive method refuge parameters for a generic Caribbean reef
- [`karpata_species`](https://cmbeese.github.io/mizerReef/reference/karpata_species.md)
  : species_params dataframe for a generic Caribbean reef
- [`karpata_model`](https://cmbeese.github.io/mizerReef/reference/karpata_model.md)
  : MizerParams object for a multispecies generic reef

### Refuge profiles

- [`aquarius_refuge`](https://cmbeese.github.io/mizerReef/reference/aquarius_refuge.md)
  : Competitive method refuge parameters for a generic Caribbean reef
