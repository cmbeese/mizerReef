# MizerReef 2.1.0

## New features

- New vignette "Combining mizerReef with mizerMR: multiple background
  resources" showing how to chain mizerReef with
  [mizerMR](https://sizespectrum.org/mizerMR/) to give the reef model several
  size-structured background resources.
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