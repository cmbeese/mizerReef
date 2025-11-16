# mizerReef 2.0.0

## Major changes
- New S4 classes: `MizerReefParams` and `MizerReefSim` for reef-specific modeling and simulation storage.
- Added slots for refuge, algae, detritus, and sponge parameters.
- Simulation now stores time-varying refuge profiles (`vulnerable`).
- Major refactor of the degradation system:
	- Removed hardcoded bleaching trajectory logic.
	- All bleaching and algae responses are now fully parameterizable via user inputs.
	- New parameters: `algae_boost`, `algae_growth_boost`, `algae_capacity_boost` for flexible post-bleaching algal dynamics.
	- Boost vectors allow compounding effects and custom duration.
	- Auto-padding and validation for boost vectors.
	- Recursive logic in `reefDegrade()` preserved, but now fully flexible.
- Improved error handling and parameter validation throughout the package.
- Documentation and examples overhauled for clarity and pkgdown compatibility.
- Lintr and roxygen2 warnings fixed.
- Vignette and article structure updated for future website integration.

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