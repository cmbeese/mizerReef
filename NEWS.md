# MizerReef 2.1.0

## New features

- New vignette "Combining mizerReef with mizerMR: multiple background
  resources" showing how to chain mizerReef with
  [mizerMR](https://sizespectrum.org/mizerMR/) to give the reef model several
  size-structured background resources.

## Bug fixes

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