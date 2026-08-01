# MizerReef 2.1.0

## New features

- New `upgradeReefParams()` function for anyone with a params object built
  with mizerReef 1.x (the version accompanying the original thesis, on this
  package's `main` branch). mizerReef 2.x is a structural rewrite: reef
  behaviour is now applied through S3 dispatch on a `mizerReef` marker class
  instead of the old direct `@rates_funcs` overrides, and reef-specific data
  (refuge, algae, detritus parameters) is laid out differently within
  `other_params`. A 1.x params object handed directly to 2.x code would
  silently run the old code path rather than erroring, so call
  `params <- upgradeReefParams(params)` once to bring it up to the current
  structure. This is a structural migration only: your tuned values (algae
  growth rate, detritus external flux, capacities, refuge method and
  parameters) are carried over as-is, not reset to package defaults --
  `upgradeReefParams()` warns you about what changed and suggests running
  `reefSteady()` afterwards if you want a fresh steady state under the
  current algae/detritus mechanism. (Objects already created with mizerReef
  2.x -- including the earlier 2.0.0 release, which stored refuge/algae/
  detritus data in dedicated S4 slots that have since been removed -- are
  upgraded automatically by `validParams()` and don't need this function.)
- New vignette "Combining mizerReef with mizerMR: multiple background
  resources" showing how to chain mizerReef with
  [mizerMR](https://sizespectrum.org/mizerMR/) to give the reef model several
  size-structured background resources.
- New `include_inverts` argument (default `FALSE`) for `plotProductivity()`
  and `plotRelativeContribution()`, so users can choose whether to include
  the "inverts" species group (previously always excluded unconditionally).
- `plotSpectraPercentChange()`/`plotlySpectraPercentChange()` are renamed to
  `plotSpectraChange()`/`plotlySpectraChange()` and gain a `use_percent`
  argument (default `TRUE`) so you can choose between a percentage change
  (e.g. `50` for a 50% increase) and the raw relative proportion (e.g.
  `0.5`), plus a `return_data` argument to get the underlying data frame
  instead of the plot.
- New `use_dummy_fish_bins` argument for `setRefuge()`/`getRefuge()`/
  `newRefuge()`/`newReefParams()` (`sigmoidal`, `binned`, and `competitive`
  refuge methods). Both settings give species-specific refuge profiles, but
  compute the underlying length-to-weight bin boundaries differently:
  - `TRUE` (default, matching this package's historical behaviour before
    this argument existed): every species' refuge length bins are converted
    to weight bins using one shared reference conversion (`a_bar`/`b_bar`),
    so all species become vulnerable/protected at the same *weight* -- but
    since species have different body shapes, the *length* at which each one
    crosses that threshold differs (reported per-species in
    `refuge_lengths`).
  - `FALSE`: each species converts the same refuge length bins to weight
    bins using its own `a`/`b`, so the weight boundaries -- and therefore
    the refuge profile in weight space -- are species-specific.

## Bug fixes

### Algae and detritus dynamics

- Algae and detritus parameters are no longer duplicated across two
  `other_params` structures. mizer's `setComponent()`/`getComponent()`
  treat `other_params$algae`/`other_params$detritus` as the single
  canonical storage location for a component's parameters, but this
  package's setters and every downstream reader/writer instead used a
  second, self-invented `other_params$algae_params`/`detritus_params`
  structure that was populated once at construction and never touched
  again -- so it silently went stale the moment `tuneUR()`/`reefSteady()`
  ran, breaking `mizer::getComponent()`/`removeComponent()` for
  algae/detritus specifically. Two related field-naming bugs in the same
  structure meant a fresh `newReefParams()` call also silently dropped
  detritus's external flux and the algae/detritus consumption matrix
  (`rho`) from the bookkeeping copy mizer itself reads. Everything now
  reads and writes the single `other_params$algae`/`other_params$detritus`
  location; the old `_params`-suffixed structures are gone. This is a pure
  storage-location consolidation -- no tuned values changed for existing
  users. `caribbean_3_model` and `caribbean_10_model` have been patched in
  place to the consolidated structure, and `upgradeReefParams()` builds it
  directly.
- Algal production no longer tracks grazer consumption. `tuneUR()`/
  `tuneUR_cc()` used to overwrite `algae_growth` with the current total
  algae consumption on every `reefSteady()` call, which silently pinned
  algae biomass constant and discarded whatever production rate was
  configured. Real algal primary production on a reef is not driven by
  grazer demand (unlike detritus, where flux-balance against consumer
  egestion is appropriate): standing algal biomass is normally low because
  grazing pressure is high, not because algae isn't growing. `algae_growth`
  is now a fixed, literature-informed constant (default `2000`
  grams/m^2/year -- see `setAlgaeParams()`'s documentation for the
  citations), and `tuneUR()`/`tuneUR_cc()` instead solve for the algae
  *biomass* that balances it, so a decrease in grazing pressure now
  increases the tuned algae biomass rather than silently reducing modelled
  production to compensate. `reefSteady()` also now leaves algae's own
  dynamics live (rather than frozen) during the fish convergence loop in
  the default case, so fish get to co-adapt to algae's real steady state
  instead of having it jumped in one shot at the end; detritus stays frozen
  throughout, since its production genuinely depends on a not-yet-tuned
  external flux. (Algae remains frozen alongside detritus when
  `new_refuge = TRUE`, matching that mode's existing behaviour of leaving
  both untouched.) `caribbean_3_model` and `caribbean_10_model` have both been
  retuned to this mechanism (`caribbean_10_model` had also never been run
  through `tuneUR()` at all before this release, so it wasn't at steady
  state under any mechanism -- the untuned original is archived at
  `inst/archive/caribbean_10_model_untuned.rda` for reference); every other
  slot, including all fish abundances, is unchanged.
- The `algae_boost`/`algae_growth_boost`/`algae_capacity_boost` feature
  (added in 2.0.0) now actually takes effect. Previously it silently did
  nothing: it tried to persist a boosted value by mutating `params` inside
  `reefDegrade()`, but `mizer::project()` treats `params` as read-only for
  the whole run, so the mutation was always discarded before the next
  timestep. The boost is now computed fresh at each timestep as a pure
  function of time, the same way `reefDegrade()` itself recomputes refuge
  density fresh rather than by mutating `params`.
- `detritus_dynamics()`/`algae_dynamics()` could produce negative, then
  `NaN`, resource biomass during a live (non-frozen) `project()` run once
  fishing shifted fish abundance far enough from the tuned steady state
  (confirmed on `caribbean_10_model` under heavy fishing effort). Fixed by
  flooring production at zero before it enters the analytic update; this
  has no effect at a `tuneUR()`/`tuneUR_cc()`-tuned steady state and only
  ever engages once a simulation has moved away from it.

### Refuge

- Fixed two bugs in `getRefuge()`'s `use_dummy_fish_bins = FALSE` branch
  (`binned` and `competitive` methods) that made it unusable: it crashed
  for any model with more than one species, and (`binned` only) refuge
  proportions were accumulated into a single vector shared across every
  species/bin combination instead of a per-species matrix, giving every
  species an identical refuge profile.
- Fixed a default-handling mismatch between `getRefuge()` and
  `reefVulnerable()` for a missing/`NULL` `use_dummy_fish_bins`: the two
  functions defaulted it in opposite directions. On the bundled
  `caribbean_3_model` -- which predates this argument -- this meant
  `reefVulnerable()` silently found no matches and left the refuge matrix
  at its all-zero initial value: **every fish was 100% vulnerable to
  predation, with the competitive refuge profile having no effect at all.**
  Both functions now default to `TRUE`, matching the package's pre-2.1.0
  behaviour; the bundled `caribbean_3_model` now stores
  `use_dummy_fish_bins = TRUE` explicitly, and `reefVulnerable()` warns if
  a mismatch like this happens again.

### Plotting functions

- `plotVulnerable()`, `plotRefugeProfile()`, `plotProductivity()`,
  `plotRelativeContribution()`, and `plotTotalBiomass()` all shared a
  species-index-misalignment bug: the logical vector used to subset
  results was recomputed *after* species had already been dropped from the
  name vector, so it no longer lined up with the original species order.
  This silently mislabelled data whenever the affected species wasn't last
  in a model's species order -- e.g. a line labelled `"herbivores"` could
  actually be plotting a different species. Separately, the first four of
  these functions also excluded any species literally named `"inverts"`
  unconditionally, with no way to opt out; `plotVulnerable()`/
  `plotRefugeProfile()` now base the exclusion on
  `species_params$refuge_user` instead of the name, and
  `plotProductivity()`/`plotRelativeContribution()` gained the
  `include_inverts` argument (see "New features"). All five functions now
  keep the species selection as a single logical vector over the original
  order throughout.
- `plotSpectraPercentChange()` (now `plotSpectraChange()`) and
  `plotProductivityRelative()`/`plotTotalBiomassRelative()`'s
  `diff_method = "percent_change"` option all promised a percentage but
  never actually multiplied by 100, so values plotted as raw fractions
  (e.g. `-0.057` instead of `-5.7`). Fixed to apply the `*100` when
  percentage output is requested (the default).
- `plotTotalAbundance()`/`plotTotalBiomass()` used a simulation's last
  *time value* directly as a *row position*, which is off by one since
  time starts at 0 -- `plotTotalAbundance()` silently plotted the
  second-to-last saved timestep, and `plotTotalBiomass()` crashed outright
  for any `mizerReefSim`. Fixed to use the actual last row, and switched
  `plotTotalBiomass()`'s species subsetting to index by name so it's
  robust to the extra algae/detritus columns mizerReef adds.
- Fixed several outright crashes and non-functional wrappers, each traced
  to leftover dead code or a copy-paste error: `plotRefugeProfile()` and
  `plotVulnerable()` crashed on `MizerSim` input (a later line silently
  undid an earlier, correct branch); `getProductivity()` -- and therefore
  `plotProductivity()`/`plotRelativeContribution()` -- never worked on a
  `MizerSim` object at all (`"object 'params' not found"`); and
  `plotlyProductivity()`/`plotlyTotalBiomassRelative()` called
  nonexistent/wrong function names instead of their intended
  `plotProductivity()`/`plotTotalBiomassRelative()` counterparts.
- `plotProductivity()`'s "no species given -> use all species" default
  (including the new `include_inverts` filtering) used `missing(species)`,
  which is `FALSE` whenever any caller in the chain forwards
  `species = species` -- true of every wrapper that calls it -- so the
  "all species" default silently never engaged through those wrappers.
  Switched to `is.null(species)`, which reflects intent correctly through
  any number of wrapper layers.

### Compatibility with other mizer extensions (e.g. mizerMR)

- The bundled `caribbean_3_model`/`caribbean_10_model` example models now
  route through mizerReef's composable rate dispatch when projected
  alongside a second stacked extension, instead of silently falling back to
  a path that assumed no other extension was present. This is what caused
  `reefSteady()` to error with "incorrect number of dimensions" when
  combining mizerReef with mizerMR.
- `getSenMort()` no longer errors when another extension changes what
  `initialNResource()` returns.
- `getBiomass()` for `mizerReefSim` objects now only folds the unstructured
  algae and detritus biomasses into the biomass table, so size-structured
  components added by other extensions (e.g. mizerMR's multiple-resource
  array) are no longer mis-reshaped into it.
- `getDegrade()` no longer errors when called on a `MizerSim` object.

### Other

- Fixed `setExtMortParams()` (and therefore `newReefParams(ext_mort_params =
  ...)`) so custom `ext_mort_params` actually works: it was previously
  broken for both a `data.frame` (crashed on first use) and a named `list`
  (failed validation outright with a misleading error) input.

# MizerReef 2.0.0

## Major changes

- New S4 class: `MizerReefParams`
	- Added slots for refuge, algae, detritus
- Refactor of the degradation system:
	- Removed hardcoded bleaching trajectory logic.
	- All bleaching and algae responses are now fully parameterizable via user inputs.
	- New parameters: `algae_boost`, `algae_growth_boost`, `algae_capacity_boost` for flexible post-bleaching algal dynamics.
	- Boost vectors allow compounding effects and custom duration.
	- Auto-padding and validation for boost vectors.
	- Recursive logic in `reefDegrade()` preserved, but now fully flexible.

- Added comprehensive vignettes and user documentation for all major workflows and features.
- All exported functions now have tests and examples to ensure reliability and reproducibility.

## Bug fixes

- Fixed double-application bug in algae growth boost during bleaching year.
- Added validation for positive numeric boost values and deg_scale structure.
- Improved type safety and maintainability in parameter handling.

## Technical notes

- All bleaching parameters are now stored in the params object for easy access.
- Compatible with existing mizerReef models (backward compatible via defaults).
- deg_scale structure unchanged: column 1 = bleaching, columns 2+ = post-bleaching.
- Recursive call pattern in `reefDegrade()` preserved (base case at t < t_bleach).
- `params@time_modified` updated on each bleaching event.

## Minor changes

- Added `progress` to Imports for progress bar support.
- Updated minimum R version to 4.1.0.
- Cleaned up example code and documentation blocks.

---
See the reference manual and pkgdown site for details on new features and usage.
