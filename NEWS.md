# MizerReef 2.1.0

## New features

- New vignette "Combining mizerReef with mizerMR: multiple background
  resources" showing how to chain mizerReef with
  [mizerMR](https://sizespectrum.org/mizerMR/) to give the reef model several
  size-structured background resources.
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