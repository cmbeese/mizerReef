# MizerReef 2.1.0

## New features

- New `upgradeReefParams()` function for anyone with a params object
  built with mizerReef 1.x (the version accompanying the original
  thesis, on this package's `main` branch). mizerReef 2.x is a
  structural rewrite -- reef behaviour is now applied through S3
  dispatch on a `mizerReef` class instead of the old direct
  `@rates_funcs` overrides, and `other_params` is laid out differently
  -- so a 1.x params object handed directly to 2.x code would silently
  run the old code path rather than erroring. Call
  `params <- upgradeReefParams(params)` once on an existing 1.x object
  to bring it up to the current structural layout. This is a
  structural migration only: your tuned values (algae growth rate,
  detritus external flux, capacities, refuge method and parameters)
  are carried over as-is, not reset to package defaults or silently
  re-tuned -- `upgradeReefParams()` warns you about what changed and
  suggests running [reefSteady()] afterwards if you want a fresh
  steady state under the current algae/detritus mechanism. (Objects
  already created with mizerReef 2.x are kept up to date automatically
  by `validParams()` and don't need this function.)
- New vignette "Combining mizerReef with mizerMR: multiple background
  resources" showing how to chain mizerReef with
  [mizerMR](https://sizespectrum.org/mizerMR/) to give the reef model several
  size-structured background resources.
- `caribbean_10_model`'s algae/detritus resource dynamics are now properly
  tuned via [tuneUR()]. Its `algae_growth` had been left at
  `newReefParams()`'s raw untuned default of `2000`, unlike
  `caribbean_3_model`'s tuned `86.87579` -- meaning algae and detritus were
  not actually at steady state for the model's calibrated abundances
  (confirmed numerically: `dB/dt` far from zero for both). The earlier,
  untuned version of this object is archived at
  `inst/archive/caribbean_10_model_untuned.rda` for reference.
- New `include_inverts` argument (default `FALSE`) for `plotProductivity()`
  and `plotRelativeContribution()`, so users can choose whether to include
  the "inverts" species group (previously always excluded unconditionally
  -- see "Bug fixes" below).
- `plotSpectraPercentChange()`/`plotlySpectraPercentChange()` are renamed to
  `plotSpectraChange()`/`plotlySpectraChange()` and gain a `use_percent`
  argument (default `TRUE`) so you can choose between a percentage change
  (e.g. `50` for a 50% increase) and the raw relative proportion (e.g.
  `0.5`), plus a `return_data` argument to get the underlying data frame
  instead of the plot. See "Bug fixes" below for why the rename happened.
- New `use_dummy_fish_bins` parameter for `setRefuge()`/`getRefuge()`/
  `newRefuge()`/`newReefParams()` (`sigmoidal`, `binned`, and `competitive`
  refuge methods). Both settings give species-specific refuge profiles, but
  compute the underlying length-to-weight bin boundaries differently:
  - `TRUE` (default, matching this package's historical behaviour before
    this parameter existed): every species' refuge length bins are
    converted to weight bins using one shared reference conversion
    (`a_bar`/`b_bar`), so all species become vulnerable/protected at the
    same *weight* -- but since species have different body shapes, the
    *length* at which each one crosses that threshold differs (reported
    per-species in `refuge_lengths`).
  - `FALSE`: each species converts the same refuge length bins to weight
    bins using its own `a`/`b`, so the weight boundaries themselves --
    and therefore the refuge profile in weight space -- are
    species-specific.

## Bug fixes

- Algae and detritus parameters are no longer duplicated across two
  `other_params` structures. `mizer::setComponent()`/`getComponent()` treat
  `other_params[[component]]` (e.g. `other_params$algae`,
  `other_params$detritus`) as the single canonical storage location for a
  component's parameters, but `setAlgaeParams()`/`setDetritusParams()` and
  every downstream reader/writer (`algae_consumption()`,
  `getAlgaeProduction()`, `getAlgaeBoost()`, `detritus_consumption()`,
  `getDetritusProduction()`, `encounter_contribution()`, `tuneUR()`/
  `tuneUR_cc()`, `scaleReefModel()`, `matchReefGrowth()`,
  `removeSpecies.mizerReef()`) instead used a second, self-invented
  `other_params$algae_params`/`detritus_params` structure. The canonical
  location was populated once at construction and never touched again, so
  it silently went stale the moment `tuneUR()`/`reefSteady()` ran --
  concretely breaking `mizer::getComponent()`/`removeComponent()`, real
  exported mizer functions, for algae/detritus specifically. Everything now
  reads and writes the single `other_params$algae`/`other_params$detritus`
  location; the `_params`-suffixed structures are gone rather than kept in
  sync. This is a pure storage-location consolidation -- no tuned values
  changed. `caribbean_3_model` and `caribbean_10_model` have been patched in
  place to the single consolidated structure, using the values each model
  was actually tuned to and running on (confirmed by each model's own
  steady-state behaviour); `upgradeReefParams()` and its tests were updated
  to build the consolidated schema.
- Algal production no longer tracks grazer consumption. `tuneUR()`/
  `tuneUR_cc()` used to overwrite `algae_growth` with the current total
  algae consumption on every `reefSteady()` call, which -- since algae
  biomass was frozen while the fish sub-model converged -- silently pinned
  algae biomass constant and discarded whatever production rate was
  configured, every time. Real algal primary production on a reef is not
  driven by grazer demand (unlike detritus, where flux-balance against
  consumer egestion is appropriate): standing algal biomass is normally low
  because grazing pressure is high, not because algae isn't growing.
  `algae_growth` is now a fixed, literature-informed constant (see
  `setAlgaeParams()`'s `algae_growth_initial` docs for the citations behind
  the default of `2000` grams/m^2/year), and `tuneUR()`/`tuneUR_cc()`
  instead solve for the algae *biomass* that balances it, so a decrease in
  grazing pressure now increases the tuned algae biomass rather than
  silently reducing modelled production to compensate. Also fixes a latent
  inconsistency where the old tuning formula used a different (feeding-level
  adjusted) consumption rate than `algae_dynamics()`/`algae_dynamics_cc()`
  themselves use, so the tuned biomass is now an exact fixed point of the
  dynamics. `reefSteady()` also no longer freezes algae's own dynamics
  during the fish convergence loop (only detritus, whose production
  genuinely depends on a not-yet-tuned external flux, is still frozen),
  since algae's dynamics now converge to a real steady state that fish need
  to be given the chance to adapt to. `caribbean_3_model` and
  `caribbean_10_model` have both been updated in place (`algae_growth` and
  algae biomass retuned to the new mechanism, plus detritus's external flux
  -- which depends on algae's raw encounter rate via egestion -- retuned to
  match; every other slot, including all fish abundances, is unchanged).
- `plotVulnerable()`, `plotRefugeProfile()`, and `plotProductivity()`
  (and therefore `plotRelativeContribution()`, which calls it) excluded any
  species literally named `"inverts"` unconditionally, with no way to opt
  out, via a hardcoded `gsub("inverts", NA, species)`. For
  `plotVulnerable()`/`plotRefugeProfile()` this is now based on
  `species_params$refuge_user` instead (species that don't use refuge have
  nothing meaningful to show on those plots -- this is what the bundled
  models' `"inverts"` group actually has `refuge_user = FALSE` for) rather
  than matching a specific name. For `plotProductivity()`/
  `plotRelativeContribution()`, a new `include_inverts` argument (default
  `FALSE`, preserving prior behaviour) lets you opt back in.

  The exclusion mechanism itself also had a real data-corruption bug in all
  four functions: the species-index vector used to subset the underlying
  result matrix was recomputed *after* the excluded species had already
  been dropped from the species-name vector, so it no longer matched the
  original species positions. This silently mislabelled data whenever the
  excluded species wasn't the last one in a model's species order (it
  happened to look correct for both bundled models purely because
  `"inverts"` is positioned last in each) -- e.g. `plotVulnerable()` could
  label a line `"herbivores"` while actually plotting a *different*
  species' vulnerability values. Fixed by keeping the selection as a single
  logical vector over the original species order throughout, rather than
  recomputing indices from an already-filtered vector.
- `plotTotalBiomass()` had the same species-index-misalignment bug
  described above, but unrelated to the `"inverts"` exclusion -- it would
  trigger for *any* explicit `species` selection that skips a species in
  the middle of a model's species order (e.g.
  `species = c("predators", "inverts")` on a 3-species model, skipping
  `"herbivores"`), silently attaching the wrong species' biomass value to
  the requested species' label. Fixed the same way, by removing the
  vestigial index-recomputation and keeping `sel_sp` as a single logical
  vector throughout.
- `plotSpectraPercentChange()`'s documentation and y-axis label promised a
  percentage (`100*(N2(w) - N1(w))/N1(w)`), but the code never actually
  multiplied by 100, so values plotted as raw fractions (e.g. `-0.057`
  instead of `-5.7`). Renamed to `plotSpectraChange()` (see "New features")
  and fixed to actually apply the `*100` when `use_percent = TRUE`
  (the default).
- `plotProductivityRelative()`/`plotTotalBiomassRelative()`'s
  `diff_method = "percent_change"` option had the identical missing-`*100`
  bug as `plotSpectraPercentChange()` above: the axis was labelled
  `"% Change in Productivity"`/`"% Change in Total Biomass"` but the value
  plotted was the raw fraction. Fixed to actually multiply by 100; the
  separate `diff_method = "rel_diff"` option (documented and implemented as
  a genuine, non-percentage relative difference) is unaffected. Also fixed
  `diff_method`'s documentation, which described the valid values as
  `` `percent.change` ``/`` `rel.diff` `` (dots) when the code has always
  only accepted `"percent_change"`/`"rel_diff"` (underscores).
- `plotTotalAbundance()`/`plotTotalBiomass()` computed `end_time <-
  max(as.numeric(dimnames(object@n)$time))` for a `MizerSim` object -- a
  *time value* -- and then used it directly as a *row position*
  (`abd[end_time, ]`). Since time starts at 0, row 1 is time 0, so this was
  off by one for every simulation: `plotTotalAbundance()` silently plotted
  the second-to-last saved timestep instead of the true last one, with no
  error or warning. `plotTotalBiomass()` had it worse -- it crashed
  outright for any `mizerReefSim` (`"arguments imply differing number of
  rows"`), because `getBiomass.mizerReefSim()` includes `algae`/`detritus`
  columns alongside the species, and the downstream species-only
  subsetting indexed that longer vector with a species-only logical vector
  instead of by name. Fixed both to use the actual last row position
  (matching the correct approach already used a few lines away in
  `plotProductivity()`'s own `MizerSim` branch), and switched
  `plotTotalBiomass()`'s species subsetting to index by name instead of by
  position so it's robust to the extra algae/detritus columns.
- `plotRefugeProfile()` crashed outright (`"no slot of name
  'species_params' for this object of class 'mizerReefSim'"`) for any
  `MizerSim` input. It correctly extracted `params <- object@params` in an
  if/else block for the `MizerSim`/`MizerParams` cases, but the very next
  line unconditionally did `params <- object` again -- dead code left over
  from a refactor that silently undid the `MizerSim` branch's correct
  extraction. Deleted the redundant line.
- `plotVulnerable()` crashed for any `MizerSim` input in the common case
  (no explicit `time_step` given): `t <- max(as.numeric(dimnames(object@n)
  $time))` correctly computed the last simulation time inside
  `if (is.null(time_step))`, but the very next line, `t <- time_step`, ran
  unconditionally and reset `t` back to `NULL`, causing
  `getVulnerable(object, time_range = t)` to fail with `"The time range
  does not contain any simulation results."` (Explicitly passing
  `time_step` happened to work, since that's the one case where the
  unconditional overwrite was actually correct.) Fixed by making the
  second assignment an `else` branch instead of a second unconditional
  statement.
- `getProductivity()` (and therefore `plotProductivity()` and
  `plotRelativeContribution()`) has never actually worked on a `MizerSim`
  object in a normal user session. Inside its `plyr::aaply()` loop over
  timesteps, it referenced a bare `params@dw`, but the local variable
  holding the params object in that branch is `sim` (`sim <- object`,
  with `sim@params` used correctly a few lines above) -- `params` was
  never defined anywhere in scope, so every call failed with `"object
  'params' not found"`. This had gone undetected because no test
  previously exercised `getProductivity()`/`plotProductivity()` on a
  `MizerSim` object, and it happened to appear to work during interactive
  development whenever an unrelated variable literally named `params`
  was sitting in scope. Fixed by referencing `sim@params@dw` directly.
- `plotlyProductivity()` called `do.call("plot2Productivity", ...)`
  instead of `do.call("plotProductivity", ...)`, and
  `plotlyTotalBiomassRelative()` called `do.call("TotalBiomassRelative",
  ...)` -- not a real function name at all (missing the `plot` prefix)
  -- instead of `do.call("plotTotalBiomassRelative", ...)`. Both were
  copy-paste errors that made the two functions entirely non-functional;
  fixed to match every other `plotly*()` wrapper's already-correct
  pattern of calling its own non-plotly counterpart.
- `plotProductivity()`'s `species` argument default handling (both for the
  pre-existing "no species given -> use all species" behaviour on the
  `MizerSim` branch, and the new `include_inverts` default filtering added
  above) used `missing(species)`, which only reflects whether an argument
  was *syntactically passed*, not whether it's meaningfully absent.
  `plot2Productivity()`/`plotProductivityRelative()`/
  `plotRelativeContribution()` all forward `species = species` to their
  inner `plotProductivity()` call unconditionally (their own `species`
  argument also defaults to `NULL`), which makes `missing(species)` FALSE
  inside `plotProductivity()` even when nobody, anywhere in the call
  chain, ever specified a species -- silently breaking the "all species"
  default (`plot2Productivity()` on two `MizerSim` objects errored
  outright with `"object 'params' not found"`) and making
  `include_inverts` a no-op through any of those three wrapper functions.
  Switched both checks to `is.null(species)`, which correctly reflects
  intent through any number of wrapper layers.
- Fixed `setExtMortParams()` (and therefore `newReefParams(ext_mort_params =
  ...)`, which calls it) so that custom `ext_mort_params` actually works.
  Previously it was completely broken for both ways a user would reasonably
  supply it: a `data.frame` passed `setExtMortParams()`'s own validation but
  was then stored as a matrix, which crashed with `"$ operator is invalid
  for atomic vectors"` the moment it was used (either immediately, inside
  `newReefParams()` itself when `include_ext_mort = TRUE`, or later in
  `reefSenMort()`); a named `list` failed `setExtMortParams()`'s own
  validation outright with a misleading "needs a column called 'nat_mort'"
  error even when `nat_mort` was provided, because converting a list with
  `as.matrix()` puts its names in the row names, not the column names.
  `setExtMortParams()` now normalises list/data-frame/matrix input to a
  plain named list before validating or storing it, matching what every
  downstream consumer already expected.
- Fixed two bugs in `getRefuge()`'s `use_dummy_fish_bins = FALSE` branch
  (`binned` and `competitive` methods) that made it unusable: `refuge_lengths`
  was built from scalar values stored in a flat named list, which crashed
  with `"length of 'dimnames' [1] not equal to array extent"` for any model
  with more than one species; and (`binned` only) the refuge proportions
  were accumulated into a single vector shared across every species/bin
  combination instead of a per-species matrix, so even without the crash,
  every species would have ended up with an identical refuge profile,
  defeating the purpose of species-specific bins.
- Fixed a default-handling mismatch between `getRefuge()` and
  `reefVulnerable()` for a missing/`NULL` `use_dummy_fish_bins`:
  `getRefuge()` defaulted it to the dummy (`TRUE`) branch via `isFALSE()`,
  but `reefVulnerable()` defaulted the same missing value to the
  species-specific (`FALSE`) branch via `isTRUE()`. On the bundled
  `caribbean_3_model` -- which predates this parameter and has no
  `use_dummy_fish_bins` field at all -- this meant `reefVulnerable()` tried
  to match `bin.id` against species-specific `"sp<i>_bin<k>"` keys that
  don't exist (its `bin.id` was built the dummy way), silently found no
  matches, and left the refuge matrix at its all-zero initial value:
  **every fish was 100% vulnerable to predation, with the competitive
  refuge profile having no effect at all.** `reefVulnerable()` now
  defaults a missing value to `TRUE`, matching `getRefuge()` and the
  package's pre-2.1.0 behaviour; the bundled `caribbean_3_model` now also
  stores `use_dummy_fish_bins = TRUE` explicitly. `reefVulnerable()` also
  now warns if `use_dummy_fish_bins = FALSE` but `bin.id` has no
  species-specific keys to match against, so a mismatch like this can't
  silently produce a no-op refuge again.
- `getSenMort()` no longer errors when another extension (such as mizerMR)
  changes what `initialNResource()` returns. Its `n_pp` default now uses
  `params@initial_n_pp` directly, which is what the internal length check
  validates against.
- `getBiomass()` for `mizerReefSim` objects now only folds the unstructured
  algae and detritus biomasses into the biomass table. Size-structured
  components added by other extensions (for example the multiple-resource
  array from mizerMR) are no longer mis-reshaped into the table, which
  previously produced a "data length differs from size of matrix" warning and
  corrupted `plotBiomass()`.
- Fixed the bundled `caribbean_3_model`/`caribbean_10_model` example models
  so that projecting them alongside a second stacked mizer extension (such as
  mizerMR) actually routes through mizerReef's composable rate dispatch,
  instead of silently falling back to a non-composable path that assumed no
  other extension was present. This is what caused `reefSteady()` to error
  with "incorrect number of dimensions" when combining mizerReef with
  mizerMR.
- `getDegrade()` no longer errors when called on a `MizerSim` object.
- The `algae_boost`/`algae_growth_boost`/`algae_capacity_boost` feature
  (added in 2.0.0) now actually takes effect. Previously it silently did
  nothing: it tried to persist a boosted value by mutating `params` inside
  `reefDegrade()`, but a single `mizer::project()` call treats `params` as
  read-only for its whole duration, so the mutation was always discarded
  before the next timestep. On top of that, the bundled example models
  stored their baseline algae growth rate under a stale field name
  (`algae_growth_initial` instead of `algae_growth`), which made the
  boost's own null-check silently fail too. The boost is now computed
  fresh at each timestep as a pure function of time by `getAlgaeBoost()`,
  applied in `getAlgaeProduction()` and `algae_dynamics_cc()` (which both
  now take a `t` argument), the same way `reefDegrade()` itself recomputes
  refuge density fresh from `t` rather than by mutating `params`.
- `detritus_dynamics()`/`algae_dynamics()` could produce negative, then
  `NaN`, resource biomass during a live (non-frozen) `project()` run.
  Their analytic ODE update is a convex combination of the current biomass
  and `production / consumption`, mathematically non-negative only if
  production is non-negative -- but detritus production includes a fixed
  `external` flux set once by `tuneUR()`/`tuneUR_cc()` at the model's
  initial abundances and never updated again, so once fishing shifted fish
  abundance enough to shrink the feces/decomposition terms, total
  production could cross zero and go negative, breaking that guarantee.
  Confirmed on `caribbean_10_model` under `effort = 1` fishing: detritus
  production went negative by year 0.3, biomass itself went negative by
  year 0.4, and the simulation eventually crashed with `NaN`. Fixed by
  flooring production at zero in both functions before it enters the
  analytic update; this has no effect at any `tuneUR()`/`tuneUR_cc()`-tuned
  steady state (production there already equals consumption, which is
  non-negative by construction) and only ever engages once a live
  simulation has moved away from it. See `detritus_dynamics()`'s and
  `algae_dynamics()`'s Details for the full mechanism.
- `caribbean_3_model`'s bundled detritus external flux (`-120.1861`, the
  same value its own model-description vignette displays) was silently
  dropped from every live simulation: it was stored under the stale field
  name `other_params$detritus_params$d_external`, left over from an old
  S4-slot-then-revert migration, while `getDetritusProduction()` reads
  `detritus_params$external` specifically -- the same class of bug already
  fixed once for algae (`algae_growth_initial` vs. `algae_growth`, above).
  Since the field was simply missing rather than present-but-wrong,
  building the production vector silently omitted it with no warning:
  the object's actual, as-shipped `dB_D/dt` was `120.1861`, not the zero
  implied by its own vignette prose and by `test-tuneUR.R`'s existing
  tests (which only ever checked a *fresh* `tuneUR(caribbean_3_model)`
  call, not the bundled object's own on-disk state). Fixed by renaming the
  field to `external`. New regression tests check the bundled objects'
  own `dB/dt` directly, not just a fresh `tuneUR()` call's output, for
  both `caribbean_3_model` and `caribbean_10_model`.
- Cleaned up stale bundled-data metadata left over from the same
  S4-slot-then-revert migration: `caribbean_3_model`/`caribbean_10_model`
  carried about a dozen dead `other_params` fields (e.g.
  `initial_algae_growth`, `initial_d_external`, `ext_decomp`,
  `sen_decomp`, `method_params`, `bin.id`, `degrade`, `refuge`,
  `refuge_lengths`, `carry_capacity`, `algae_capacity`,
  `detritus_capacity`) at the top level of `other_params` -- remnants of
  an old flat-storage convention, superseded by the current nested
  `algae_params`/`detritus_params`/`refuge_params` sub-lists but never
  removed. Confirmed unread by any current code path and removed from
  both bundled models, which now match a freshly-built model's
  `other_params` structure exactly. `caribbean_10_model` also carried a
  dead duplicate `d_external` field alongside the correct `external`
  field (from before it was `tuneUR()`'d); removed. Two vignettes
  (`karpata_model-description.Rmd`, `caribbean_3_model-description.Rmd`)
  displayed prose numbers read from these stale top-level fields instead
  of the correct nested ones; repointed to the correct fields (the
  displayed values are unchanged, since both copies held identical data).
- `newReefParams()` populated `mizer::setComponent()`'s bookkeeping copy
  of the algae/detritus components (`other_params$algae`,
  `other_params$detritus`) from field names that no setter in the current
  codebase writes to -- `algae_growth_initial` and `d_external` for the
  default `use_UR_cc = FALSE` path (the same stale names as the two bugs
  above), and a set of flat top-level fields (`algae_capacity`,
  `initial_algae_growth`, `sen_decomp`, `ext_decomp`, `detritus_capacity`,
  `external`) for `use_UR_cc = TRUE`, left over from `main`'s old
  flat-storage convention. A freshly-built model's bookkeeping copy was
  therefore silently `NULL` on both paths. The two bundled models only
  ever displayed correct numbers from this copy because their `.rda`
  files were built with older, still-self-consistent code and the copy
  is set once at construction and never refreshed -- a frozen historical
  accident, not evidence the current code path works. No simulation
  dynamics are affected (`algae_dynamics()`/`detritus_dynamics()` and
  their `_cc` variants all read `other_params$algae_params`/
  `detritus_params` directly, never this bookkeeping copy), but display
  code that reads it -- as `caribbean_3_model-description.Rmd` and
  `karpata_model-description.Rmd` both did, for algal/detrital
  production prose -- would silently break on a freshly-rebuilt model.
  Fixed `newReefParams()` to read the current field names on both paths,
  and repointed both vignettes to read the authoritative nested fields
  directly rather than this bookkeeping copy (displayed numbers are
  unchanged).
- `newReefParams()` stored the algae/detritus consumption matrix under
  `other_params$algae_params$rho_algae`/
  `other_params$detritus_params$rho_detritus`, but every consumer
  (`algae_components.R`, `detritus_components.R`, `reef-components.R`,
  `reef-helpers.R`) reads and writes the bare field name `rho`. This only
  worked because R's `$` operator silently does partial name-matching
  (`rho` uniquely prefix-matches `rho_algae`) -- the same migration-gap
  bug class as the `d_external`/`algae_growth_initial` fixes above, just
  one that hadn't broken output yet. Fixed `newReefParams()` to write
  `rho` directly, matching every reader, and patched both bundled models
  (`caribbean_3_model`, `caribbean_10_model`) the same way. New
  regression tests guard the exact field name on both bundled models and
  on a freshly-built model.
- Removed a few more remnants of the same stale-bundled-data audit:
  `caribbean_10_model` carried a dead duplicate `new_refuge` field nested
  inside `refuge_params` (the live flag lives at the top level of
  `other_params`) and an extra `capacity` field in its `algae`
  bookkeeping copy that `caribbean_3_model` lacked -- both confirmed
  unread by any current code path and removed so the two bundled models'
  `other_params` structures now match exactly.
- Removed a few lines of dead code in `getDegrade()`/`getVulnerable()`
  (`R/reef-rates.R`) that read `other_params$dt` and the stale top-level
  `other_params$method_params` -- neither value was ever used after being
  assigned.

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