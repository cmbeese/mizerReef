# Changelog

## MizerReef 2.0.3

### New features

- The reports mizerReef makes while building or changing a model go
  through mizer 3.3’s
  [`signal_info()`](https://sizespectrum.org/mizer/reference/signal_info.html)/[`with_info_level()`](https://sizespectrum.org/mizer/reference/with_info_level.html)
  mechanism instead of bare
  [`message()`](https://rdrr.io/r/base/message.html)/[`warning()`](https://rdrr.io/r/base/warning.html)
  calls, so they can be turned down.
  [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md),
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md),
  [`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md),
  [`setDetritusParams()`](https://cmbeese.github.io/mizerReef/reference/setDetritusParams.md),
  [`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md),
  [`newRefuge()`](https://cmbeese.github.io/mizerReef/reference/newRefuge.md),
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)
  and
  [`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
  gain an `info_level` argument, defaulting to
  [`mizer::default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.html):

  ``` r

  params <- newReefParams(..., info_level = 1)   # only what went differently
  params <- reefSteady(params, info_level = 0)   # silence
  ```

  The wording is unchanged, and every report is still a warning at the
  default level, so scripts matching on these messages keep working –
  except
  [`newRefuge()`](https://cmbeese.github.io/mizerReef/reference/newRefuge.md)’s
  “no new method given” report: it was a single multi-line
  [`warning()`](https://rdrr.io/r/base/warning.html) string literal, so
  its text included a literal newline plus the source file’s own
  indentation; the [`paste()`](https://rdrr.io/r/base/paste.html)-built
  replacement collapses that to a single space. No test matched the old
  text exactly, so nothing in this package’s own suite broke, but a
  caller doing an exact-string match against it would need updating.

  [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
  and
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  install one handler for the whole call, so the reports raised by the
  setters and by
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
  are collected and given together. `info_level` is also forwarded to
  those inner calls explicitly. Nested
  [`with_info_level()`](https://sizespectrum.org/mizer/reference/with_info_level.html)
  calls do work without forwarding it in the common case – but only
  because an inner call’s own resolved `info_level` (read from
  [`mizer::default_info_level()`](https://sizespectrum.org/mizer/reference/default_info_level.html))
  happens to agree with the outer one, which is true whenever both are
  left at the default. If a global `options(mizer_info_level = 0)`
  differs from what the outer call was explicitly given, an unforwarded
  inner call resolves its own default independently, gets `0`, and
  [`with_info_level()`](https://sizespectrum.org/mizer/reference/with_info_level.html)‘s
  documented “silence is the exception” rule takes over: it
  unconditionally muffles that call’s reports regardless of what the
  outer, explicit `info_level` asked for. Confirmed directly: under a
  global `mizer_info_level = 0`, `reefSteady(params, info_level = 3)`
  silently lost
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)’s
  reports even though the caller explicitly asked for them; explicit
  forwarding fixes it. There is no argument-collision risk in forwarding
  it – `info_level` is one of these functions’ own named formals, so it
  can never also be present in `...` for R to match twice, contrary to
  what was assumed here previously.

- [`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md)
  reports the `r_max` -\> `R_max` rename with
  [`signal_info()`](https://sizespectrum.org/mizer/reference/signal_info.html),
  matching mizer’s own
  [`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.html),
  which this part of the function reproduces. Note that, also matching
  upstream, this branch is unreachable in practice:
  [`validParams()`](https://sizespectrum.org/mizer/reference/validParams.html),
  called earlier in the same function, already silently renames `r_max`
  to `R_max` on its own before this check runs.

- [`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)’s
  report for missing `a`/`b` length-weight parameters now identifies
  which column was actually defaulted (`var = "a"` or `"b"`) instead of
  always reporting `"a"`, so a caller suppressing this specific report
  by name via
  [`with_info_level()`](https://sizespectrum.org/mizer/reference/with_info_level.html)’s
  `except` argument targets the right one.

### Not converted, deliberately

Ten reports stay as plain
[`message()`](https://rdrr.io/r/base/message.html)/[`warning()`](https://rdrr.io/r/base/warning.html)
calls:

- The `is.nan(consumption)` guards in
  [`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md),
  [`algae_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics_cc.md),
  [`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md)
  and
  [`detritus_dynamics_cc()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics_cc.md),
  and the refuge-misconfiguration warning in
  [`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md).
  These run inside
  [`project()`](https://sizespectrum.org/mizer/reference/project.html)’s
  per-timestep loop, where there is no wrapped entry point to collect
  them and end-of-call delivery would be the wrong semantics. They are
  also genuine fault alarms rather than a choice mizer made, which is
  what
  [`signal_info()`](https://sizespectrum.org/mizer/reference/signal_info.html)
  is documented for.
- The four migration reports in
  [`upgradeReefParams()`](https://cmbeese.github.io/mizerReef/reference/upgradeReefParams.md),
  which fire once during a deliberate one-off migration.
- The warning in
  [`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md):
  a plotting function has no `info_level` convention.

## MizerReef 2.0.2

Upgraded to mizer 3.3. `DESCRIPTION` now requires `mizer (>= 3.3.0)` and
`mizerExperimental (>= 3.3.0)`.

### Bug fixes

- [`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)‘s
  `a`/`b` default-fill (from `a_bar`/`b_bar`, for species missing their
  own length-weight parameters) is kept as a direct `@species_params`
  slot write, unlike the rest of this function and unlike every other
  conversion in this release. Routing it through
  `species_params(params, recalculate = FALSE) <-` – even with
  `recalculate = FALSE` – runs mizer’s built-in
  `a`/`b`-vs-`l_max`/`w_max` consistency check, which silently
  overwrites `l_max` to match whenever `a * l_max^b` disagrees with the
  model’s real `w_max`. A generic `a_bar`/`b_bar` fallback disagreeing
  with a specific species’ real `l_max`/`w_max` is the normal case, not
  an edge case – it’s exactly the situation this fallback exists for.
  Confirmed directly: filling missing `a`/`b` with `a_bar = 0.04`,
  `b_bar = 3.2` on `caribbean_3_model` silently moved `l_max` from 50 to
  33.8 (32%) when routed through the setter, and the same corruption
  resurfaces at the *next* full-table `species_params<-()` write
  anywhere downstream
  (e.g. [`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)/
  [`setDetritusParams()`](https://cmbeese.github.io/mizerReef/reference/setDetritusParams.md),
  both called right after
  [`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)
  inside
  [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md))
  for as long as the mismatched `a`/`b` remain stored. `a`/`b` filled
  this way are therefore not recorded as given – a subsequent
  recalculation could revert them to `NA` – but that is far more benign
  than silently corrupting `l_max`.

- [`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md)’s
  new `biomass` and `per_log_size` arguments were inserted between
  `power` and `use_percent`, silently reinterpreting old positional
  calls: `power` used to be required (so commonly supplied
  positionally), and any caller also supplying `use_percent`
  positionally (previously the 5th argument, now the 7th) had that value
  silently reinterpreted as `biomass` instead, with `use_percent`
  quietly falling back to its default. Reordered so `use_percent` keeps
  its original position and the two new arguments come after it.

- `plotSpectraChange(..., size_axis = "l")` joined `object1`‘s and
  `object2`’s spectra on their length columns, but length is derived
  per-object from that object’s own `a`/`b` – so whenever the two
  objects being compared don’t share identical `a`/`b` (the normal case:
  comparing two different objects is this function’s whole purpose), the
  two length columns disagree at the same underlying size, and the join
  silently dropped those rows as `NA` instead of erroring. Confirmed:
  comparing two models differing only in a 10% change to species `a`
  produced 300 `NA` rows out of 376. The join now always happens on `w`
  – the model’s shared discretisation grid – with `object1`’s own length
  values attached afterwards for display when a length axis is
  requested, so no rows are dropped regardless of how the two objects’
  `a`/`b` differ.

- [`matchReefGrowth()`](https://cmbeese.github.io/mizerReef/reference/matchReefGrowth.md)
  records the parameters it scales, so a later recalculation can no
  longer undo the growth match. It scales
  `search_vol`/`intake_max`/`metab` by hand and scales
  `gamma`/`h`/`ks`/`k` to match, but wrote those scalars into the
  `@species_params` slot. That left
  [`species_params()`](https://sizespectrum.org/mizer/reference/species_params.html)
  disagreeing with
  [`given_species_params()`](https://sizespectrum.org/mizer/reference/species_params.html),
  and the next call that triggered a recalculation restored the unscaled
  values: on a model built this way, a no-op
  `species_params(p) <- species_params(p)` moved the metabolic rate by
  88% and the search volume by 95%. The scalars now go in through
  `species_params(params, recalculate = FALSE) <-`, which is what
  mizer’s own
  [`matchGrowth()`](https://sizespectrum.org/mizer/reference/matchGrowth.html)
  does and for the same reason. A freshly built reef model is now
  unchanged by a recalculation.

  **The bundled `caribbean_3_model` and `caribbean_10_model` were built
  with the old code and still carry the inconsistency** (`given` `ks` of
  0.15 vs a `used` 0.0797 in `caribbean_3_model`). They are unaffected
  in ordinary use, but any call that recalculates the species parameters
  will shift their rates. They need rebuilding from `inst/scripts/` to
  be fully consistent.

- Every species parameter mizerReef sets now goes in through
  `species_params(params, recalculate = FALSE) <-` rather than by
  writing the `params@species_params` slot: `rho_algae`/`rho_detritus`
  in
  [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md),
  [`rescale_algae()`](https://cmbeese.github.io/mizerReef/reference/rescale_algae.md),
  [`rescale_detritus()`](https://cmbeese.github.io/mizerReef/reference/rescale_detritus.md)
  and
  [`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md);
  the `interaction_<resource>` columns in
  [`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)/[`setDetritusParams()`](https://cmbeese.github.io/mizerReef/reference/setDetritusParams.md);
  `a`/`b`, `refuge_user`, `blocked_pred` and `satiation` in
  [`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md);
  and the `bad_pred` -\> `blocked_pred` rename in
  [`upgradeReefParams()`](https://cmbeese.github.io/mizerReef/reference/upgradeReefParams.md).
  The values are unchanged – `recalculate = FALSE` rebuilds nothing,
  exactly as the slot write did – but the table is now validated and the
  values are recorded as given, so mizer will not treat them as its own
  defaults and recalculate over them.

  One slot write is deliberately kept: the `constant_reproduction` flag
  in
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md),
  which is set and removed within the one call exactly as mizer’s own
  `tune_steady_project()` does it.

- [`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md)’s
  block that reproduces mizer’s
  [`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.html)
  also wrote its scaled `R_max` and `gamma` directly into the
  `@species_params` slot, the same bug as above, inherited from mizer’s
  own
  [`scaleModel()`](https://sizespectrum.org/mizer/reference/scaleModel.html)
  (confirmed against mizer 3.3.1: `scaleModel(NS_params, factor = 2)`
  leaves `given_species_params()$R_max` untouched while
  [`species_params()`](https://sizespectrum.org/mizer/reference/species_params.html)
  doubles, and a later no-op recalculation silently reverts it –
  reported upstream as
  [sizespectrum/mizer#599](https://github.com/sizespectrum/mizer/issues/599)).
  [`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md)
  now routes these through
  `species_params(params, recalculate = FALSE) <-` too, rather than
  waiting on the upstream fix, so
  [`calibrateReefBiomass()`](https://cmbeese.github.io/mizerReef/reference/calibrateReefBiomass.md)
  and
  [`calibrateReefNumber()`](https://cmbeese.github.io/mizerReef/reference/calibrateReefNumber.md)
  (which both call
  [`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md))
  no longer risk this drift.

- [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)’s
  `...` was spliced into its
  [`mizer::findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.html)/
  [`mizer::projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.html)
  calls alongside arguments those calls already hardcode by name
  (`distance_func`, `t_check`, `t_save`, `distance_tol`,
  `require_steady`, `solver`). Passing any of those names through `...`
  – e.g. `reefSteady(params, distance_tol = 0.05)`, or the same via
  [`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)/[`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.html)
  – errored with `formal argument matched by multiple actual arguments`
  instead of overriding the value. The defaults are now merged with
  `...` via [`modifyList()`](https://rdrr.io/r/utils/modifyList.html)
  (user values win) instead of being passed alongside it, so any of
  these mizer 3.3 arguments can be overridden as documented.

### Breaking changes

- [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  is no longer written over
  [`mizer::steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
  in mizer’s namespace. It is registered as the
  [`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
  **and**
  [`tuneSteadyState()`](https://sizespectrum.org/mizer/reference/tuneSteadyState.html)
  method for `mizerReef` objects instead, so `steady(reef_params)` and
  `tuneSteadyState(reef_params)` both do the reef-aware thing while
  every non-reef model in the session is left alone.

  The old `utils::assignInNamespace("steady", reefSteady, ns = "mizer")`
  replaced mizer’s S3 *generic*, not just its `MizerParams` method.
  Under mizer 3.3 that had three consequences:
  [`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
  on any non-reef model errored with `argument is of length zero`
  (reefSteady() reads `params@other_params$new_refuge`, which is `NULL`
  there); no other extension package’s
  [`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
  method could ever dispatch; and mizer 3.3’s new
  [`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
  arguments (`t_save`, `amplitude_tol`, `amp_rel_tol`,
  `extinction_threshold`, `info_level`, `method`) were silently
  swallowed. Note that a session that loaded an older mizerReef keeps
  the patched
  [`mizer::steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
  until R is restarted.

- [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  passes `...` on to mizer’s steady-state finder rather than to
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md),
  which ignored it. This is what lets `effort`, `method` and
  `info_level` reach the projection. A misspelled argument that used to
  be silently discarded is now an error.

- [`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md)
  labels its y axis after the quantity actually plotted, so the default
  is now `"% Change in Biomass density"` rather than
  `"% Change in Biomass"`.

- [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
  and
  [`upgradeReefParams()`](https://cmbeese.github.io/mizerReef/reference/upgradeReefParams.md)
  record only mizerReef in `@extensions`, via
  [`mizer::recordExtension()`](https://sizespectrum.org/mizer/reference/recordExtension.html),
  instead of copying the whole session extension registry with
  [`getRegisteredExtensions()`](https://sizespectrum.org/mizer/reference/getRegisteredExtensions.html).
  An extension can be loaded without having been applied to a particular
  model, so the old code made the model claim it: with mizerMR merely
  loaded,
  [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
  returned an object recording mizerMR and promoted to S4 class
  `mizerMR`, with no `MR` component for mizerMR’s methods to read.
  Chaining is unaffected – each extension records itself as it is
  applied, and the bundled `caribbean_3_model` has only ever recorded
  mizerReef – so
  [`setMultipleResources()`](https://sizespectrum.org/mizerMR/reference/setMultipleResources.html)
  still produces a properly chained `mizerMR`/`mizerReef` object.

### New features

- [`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md)
  gained the `biomass` and `per_log_size` arguments that mizer 3.3
  introduced in place of `power`, and `power` is now optional,
  defaulting to the biomass density as in
  [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.html).
  Passing `biomass` through `...` alongside `power` previously reached
  [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.html)
  as a contradictory pair, which mizer 3.3 rejects.

### Internal

- [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  calls
  [`mizer::findSteadyState()`](https://sizespectrum.org/mizer/reference/findSteadyState.html)
  (or
  [`mizer::projectUntilSettled()`](https://sizespectrum.org/mizer/reference/projectUntilSettled.html)
  when `return_sim = TRUE`) instead of the superseded
  [`mizer::projectToSteady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html).
  It asks for the old stopping rule explicitly
  (`require_steady = FALSE`), exactly as mizer’s own
  [`steady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
  and
  [`projectToSteady()`](https://sizespectrum.org/mizer/reference/superseded_steady.html)
  wrappers do, so convergence behaviour is unchanged:
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  still stops on its documented `tol` criterion and does not
  additionally require the biomass drift to settle.

- Dropped the `splus2R` dependency. It was used for a single
  [`splus2R::is.number()`](https://rdrr.io/pkg/splus2R/man/is.number.html)
  call in
  [`getSenMort()`](https://cmbeese.github.io/mizerReef/reference/getSenMort.md),
  where every other predicate in the same
  [`assert_that()`](https://rdrr.io/pkg/assertthat/man/assert_that.html)
  came from `assertthat`, which the package already imports and whose
  [`is.number()`](https://rdrr.io/pkg/assertthat/man/scalar.html) makes
  the same check.

- Documentation and vignettes updated for mizer 3.3: `calibrateYield()`
  has been removed from mizer, and mizer’s extension articles were
  renamed to `guide-use-extension-packages` and
  `guide-create-extension-package`. Example code now uses
  `biomass`/`per_log_size` in place of `power`. This needs mizerMR
  0.3.1.2 or later in
  [`vignette("using-multiple-resources")`](https://cmbeese.github.io/mizerReef/articles/using-multiple-resources.md):
  0.3.1.1’s
  [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.html)
  method still had the pre-3.3 signature and passed its own `power` down
  to [`NextMethod()`](https://rdrr.io/r/base/UseMethod.html), so the
  flags reached mizer as a contradictory pair on a mizerReef+mizerMR
  model. `Suggests` now requires `mizerMR (>= 0.3.1.2)`.

- [`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md)
  names its y axis correctly on a model whose
  [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.html)
  method renames the value column generically, as mizerMR’s does (it
  calls it `value`, so the label came out as `"% Change in value"`). The
  plotted quantity is then derived from `power`/`biomass`/`per_log_size`
  using mizer’s own rule instead.

- [`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md)
  and
  [`plotlySpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md)
  follow
  [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.html)
  onto a length axis: `size_axis = "l"` used to error in the join,
  because
  [`plotSpectra()`](https://sizespectrum.org/mizer/reference/plotSpectra.html)
  returns an `l` column rather than a `w` one.

## MizerReef 2.0.1

### Bug fixes

- [`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md)
  (called by
  [`calibrateReefBiomass()`](https://cmbeese.github.io/mizerReef/reference/calibrateReefBiomass.md)
  and
  [`calibrateReefNumber()`](https://cmbeese.github.io/mizerReef/reference/calibrateReefNumber.md))
  no longer rescales `algae_growth`. The 2.0.0 algae redesign (see
  “Algae and detritus dynamics” under 2.0.0 below) fixed `algae_growth`
  as a constant, literature-informed production rate that
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
  solve algae’s biomass around, but
  [`scaleReefModel()`](https://cmbeese.github.io/mizerReef/reference/scaleReefModel.md)
  was never updated to match – it still applied the same
  fish-abundance-correction factor to `algae_growth` that it correctly
  applies to fish abundance and to the algae/detritus encounter
  coefficients (`rho`/`rho_algae`/`rho_detritus`). Since nothing
  downstream ever resets `algae_growth` afterwards, a single
  [`calibrateReefBiomass()`](https://cmbeese.github.io/mizerReef/reference/calibrateReefBiomass.md)
  call against a starting point far from the observed target could
  silently collapse it by many orders of magnitude, corrupting algae’s
  production for the rest of the calibration pipeline.
  `rho`/`rho_algae`/`rho_detritus` continue to rescale as before – that
  part keeps algae/detritus consumption pressure invariant under the
  abundance rescale, the same role `search_vol`/`gamma` play for
  fish-fish encounters, and was already correct.
- Detritus and algae colors are now correctly set in
  [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
  to not be overwritten by
  [`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html).
  That fix’s own re-assignment, plus the equivalent assignments in
  [`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)/[`setDetritusParams()`](https://cmbeese.github.io/mizerReef/reference/setDetritusParams.md),
  unconditionally wrote `algae_colour`/`detritus_colour` into
  `@linecolour` even when `NULL` – a documented, valid “leave it as set”
  input for both – which throws `replacement has length zero` on a plain
  character vector. All three now skip the assignment when the colour is
  `NULL`, matching how every other optional parameter in these functions
  already handles `NULL`.

## MizerReef 2.0.0

This is the first public release of mizerReef beyond 1.0.1, the version
that accompanied the original thesis. Everything below describes how
2.0.0 differs from that 1.0.1 release, which is what any existing user
will actually have installed – there is no released version in between
for these notes to assume you’ve already seen.

### New features

- New
  [`upgradeReefParams()`](https://cmbeese.github.io/mizerReef/reference/upgradeReefParams.md)
  function for anyone with a params object built with mizerReef 1.x (the
  version accompanying the original thesis, on this package’s `main`
  branch). mizerReef 2.x is a structural rewrite: reef behaviour is now
  applied through S3 dispatch on a `mizerReef` marker class instead of
  the old direct `@rates_funcs` overrides, and reef-specific data
  (refuge, algae, detritus parameters) is laid out differently within
  `other_params`. A 1.x params object handed directly to 2.x code would
  silently run the old code path rather than erroring, so call
  `params <- upgradeReefParams(params)` once to bring it up to the
  current structure. This is a structural migration only: your tuned
  values (algae growth rate, detritus external flux, capacities, refuge
  method and parameters) are carried over as-is, not reset to package
  defaults –
  [`upgradeReefParams()`](https://cmbeese.github.io/mizerReef/reference/upgradeReefParams.md)
  warns you about what changed and suggests running
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  afterwards if you want a fresh steady state under the current
  algae/detritus mechanism. (This function is only needed for 1.x
  objects; any object already created with this 2.0.0 release stays
  current automatically via
  [`validParams()`](https://sizespectrum.org/mizer/reference/validParams.html).)
- New vignette “Combining mizerReef with mizerMR: multiple background
  resources” showing how to chain mizerReef with
  [mizerMR](https://sizespectrum.org/mizerMR/) to give the reef model
  several size-structured background resources.
- New `include_inverts` argument (default `FALSE`) for
  [`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)
  and
  [`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md),
  so users can choose whether to include the “inverts” species group
  (previously always excluded unconditionally).
- `plotSpectraPercentChange()`/`plotlySpectraPercentChange()` are
  renamed to
  [`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md)/[`plotlySpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md)
  and gain a `use_percent` argument (default `TRUE`) so you can choose
  between a percentage change (e.g. `50` for a 50% increase) and the raw
  relative proportion (e.g. `0.5`), plus a `return_data` argument to get
  the underlying data frame instead of the plot.
- New `use_dummy_fish_bins` argument for
  [`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)/[`getRefuge()`](https://cmbeese.github.io/mizerReef/reference/getRefuge.md)/
  [`newRefuge()`](https://cmbeese.github.io/mizerReef/reference/newRefuge.md)/[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
  (`sigmoidal`, `binned`, and `competitive` refuge methods). Both
  settings give species-specific refuge profiles, but compute the
  underlying length-to-weight bin boundaries differently:
  - `TRUE` (default, matching this package’s historical behaviour before
    this argument existed): every species’ refuge length bins are
    converted to weight bins using one shared reference conversion
    (`a_bar`/`b_bar`), so all species become vulnerable/protected at the
    same *weight* – but since species have different body shapes, the
    *length* at which each one crosses that threshold differs (reported
    per-species in `refuge_lengths`).
  - `FALSE`: each species converts the same refuge length bins to weight
    bins using its own `a`/`b`, so the weight boundaries – and therefore
    the refuge profile in weight space – are species-specific.
- [`setRefuge()`](https://cmbeese.github.io/mizerReef/reference/setRefuge.md)/[`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)’s
  `satiation` argument (which controls which species get a
  satiation-limited feeding response, versus unlimited intake) is now
  optional. 1.x required every model to specify `satiation` explicitly
  and errored otherwise; it’s now defaulted automatically – `TRUE` only
  for species whose diet is detritus without also grazing algae or
  eating other species (pure detritivores), `FALSE` for everyone else,
  including herbivores that graze algae, carnivores, and species that
  eat both algae and detritus. This matches the design intent that
  satiation-limited consumption applies to detritivory specifically, not
  herbivory in general – a species that grazes algae behaves like an
  unregulated herbivore even if it also scavenges some detritus. A
  warning is issued whenever the default is used, so it’s always visible
  when a model is relying on it rather than an explicit choice.
- Habitat degradation (bleaching) trajectories are now fully
  user-parameterizable. 1.x had a single bleaching-severity trajectory
  hardcoded directly into the rate-calculation code, with no way to
  configure it without editing the package source.
  [`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
  now takes a `deg_scale` trajectory matrix/data.frame and
  `algae_boost`/ `algae_growth_boost`/`algae_capacity_boost` vectors as
  ordinary arguments, so post-bleaching algae/refuge dynamics can be
  customised per model without touching package code.

### Bug fixes

#### Algae and detritus dynamics

- Algae and detritus parameters are no longer duplicated across two
  `other_params` structures. mizer’s
  [`setComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)/[`getComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)
  treat `other_params$algae`/`other_params$detritus` as the single
  canonical storage location for a component’s parameters, but this
  package’s setters and every downstream reader/writer instead used a
  second, self-invented `other_params$algae_params`/`detritus_params`
  structure that was populated once at construction and never touched
  again – so it silently went stale the moment
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  ran, breaking
  [`mizer::getComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)/[`removeComponent()`](https://sizespectrum.org/mizer/reference/setComponent.html)
  for algae/detritus specifically. Two related field-naming bugs in the
  same structure meant a fresh
  [`newReefParams()`](https://cmbeese.github.io/mizerReef/reference/newReefParams.md)
  call also silently dropped detritus’s external flux and the
  algae/detritus consumption matrix (`rho`) from the bookkeeping copy
  mizer itself reads. Everything now reads and writes the single
  `other_params$algae`/`other_params$detritus` location; the old
  `_params`-suffixed structures are gone. This is a pure
  storage-location consolidation – no tuned values changed for existing
  users. `caribbean_3_model` and `caribbean_10_model` have been patched
  in place to the consolidated structure, and
  [`upgradeReefParams()`](https://cmbeese.github.io/mizerReef/reference/upgradeReefParams.md)
  builds it directly.
- Algal production no longer tracks grazer consumption.
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/
  [`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
  used to overwrite `algae_growth` with the current total algae
  consumption on every
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  call, which silently pinned algae biomass constant and discarded
  whatever production rate was configured. Real algal primary production
  on a reef is not driven by grazer demand (unlike detritus, where
  flux-balance against consumer egestion is appropriate): standing algal
  biomass is normally low because grazing pressure is high, not because
  algae isn’t growing. `algae_growth` is now a fixed,
  literature-informed constant (default `2000` grams/m^2/year – see
  [`setAlgaeParams()`](https://cmbeese.github.io/mizerReef/reference/setAlgaeParams.md)’s
  documentation for the citations), and
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)
  instead solve for the algae *biomass* that balances it, so a decrease
  in grazing pressure now increases the tuned algae biomass rather than
  silently reducing modelled production to compensate.
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  also now leaves algae’s own dynamics live (rather than frozen) during
  the fish convergence loop in the default case, so fish get to co-adapt
  to algae’s real steady state instead of having it jumped in one shot
  at the end; detritus stays frozen throughout, since its production
  genuinely depends on a not-yet-tuned external flux. (Algae remains
  frozen alongside detritus when `new_refuge = TRUE`, matching that
  mode’s existing behaviour of leaving both untouched.)
  `caribbean_3_model` and `caribbean_10_model` have both been retuned to
  this mechanism (`caribbean_10_model` had also never been run through
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)
  at all before this release, so it wasn’t at steady state under any
  mechanism – the untuned original is archived at
  `inst/archive/caribbean_10_model_untuned.rda` for reference); every
  other slot, including all fish abundances, is unchanged.
- The `algae_boost`/`algae_growth_boost`/`algae_capacity_boost` feature
  now actually takes effect. Previously it silently did nothing: it
  tried to persist a boosted value by mutating `params` inside
  [`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md),
  but
  [`mizer::project()`](https://sizespectrum.org/mizer/reference/project.html)
  treats `params` as read-only for the whole run, so the mutation was
  always discarded before the next timestep. The boost is now computed
  fresh at each timestep as a pure function of time, the same way
  [`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md)
  itself recomputes refuge density fresh rather than by mutating
  `params`.
- [`detritus_dynamics()`](https://cmbeese.github.io/mizerReef/reference/detritus_dynamics.md)/[`algae_dynamics()`](https://cmbeese.github.io/mizerReef/reference/algae_dynamics.md)
  could produce negative, then `NaN`, resource biomass during a live
  (non-frozen)
  [`project()`](https://sizespectrum.org/mizer/reference/project.html)
  run once fishing shifted fish abundance far enough from the tuned
  steady state (confirmed on `caribbean_10_model` under heavy fishing
  effort). Fixed by flooring production at zero before it enters the
  analytic update; this has no effect at a
  [`tuneUR()`](https://cmbeese.github.io/mizerReef/reference/tuneUR.md)/[`tuneUR_cc()`](https://cmbeese.github.io/mizerReef/reference/tuneUR_cc.md)-tuned
  steady state and only ever engages once a simulation has moved away
  from it.
- [`getDetritusProduction()`](https://cmbeese.github.io/mizerReef/reference/getDetritusProduction.md)’s
  documented formula for the “decomp” term applied a single, unnamed
  `prop_decomp` weight to only the external-mortality contribution and
  none to the senescence-mortality contribution, and its “feces” term’s
  integral omitted the abundance factor `N_i(w)` entirely. Neither
  matched the actual code, which weights the senescence and
  external-mortality contributions by `sen_decomp` and `ext_decomp`
  respectively (as
  [`setDetritusParams()`](https://cmbeese.github.io/mizerReef/reference/setDetritusParams.md)’s
  own documentation already correctly describes) and does integrate over
  abundance. Documentation fixed to match actual behaviour; no code
  change.

#### Refuge

- Fixed two bugs in
  [`getRefuge()`](https://cmbeese.github.io/mizerReef/reference/getRefuge.md)’s
  `use_dummy_fish_bins = FALSE` branch (`binned` and `competitive`
  methods) that made it unusable: it crashed for any model with more
  than one species, and (`binned` only) refuge proportions were
  accumulated into a single vector shared across every species/bin
  combination instead of a per-species matrix, giving every species an
  identical refuge profile.
- Fixed a default-handling mismatch between
  [`getRefuge()`](https://cmbeese.github.io/mizerReef/reference/getRefuge.md)
  and
  [`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
  for a missing/`NULL` `use_dummy_fish_bins`: the two functions
  defaulted it in opposite directions. On the bundled
  `caribbean_3_model` – which predates this argument – this meant
  [`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
  silently found no matches and left the refuge matrix at its all-zero
  initial value: **every fish was 100% vulnerable to predation, with the
  competitive refuge profile having no effect at all.** Both functions
  now default to `TRUE`, matching the package’s intended original
  behaviour; the bundled `caribbean_3_model` now stores
  `use_dummy_fish_bins = TRUE` explicitly, and
  [`reefVulnerable()`](https://cmbeese.github.io/mizerReef/reference/reefVulnerable.md)
  warns if a mismatch like this happens again.

#### Plotting functions

- [`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md),
  [`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md),
  [`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
  [`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md),
  and
  [`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)
  all shared a species-index-misalignment bug: the logical vector used
  to subset results was recomputed *after* species had already been
  dropped from the name vector, so it no longer lined up with the
  original species order. This silently mislabelled data whenever the
  affected species wasn’t last in a model’s species order – e.g. a line
  labelled `"herbivores"` could actually be plotting a different
  species. Separately, the first four of these functions also excluded
  any species literally named `"inverts"` unconditionally, with no way
  to opt out;
  [`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)/
  [`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md)
  now base the exclusion on `species_params$refuge_user` instead of the
  name, and
  [`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)/[`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md)
  gained the `include_inverts` argument (see “New features”). All five
  functions now keep the species selection as a single logical vector
  over the original order throughout.
- `plotSpectraPercentChange()` (now
  [`plotSpectraChange()`](https://cmbeese.github.io/mizerReef/reference/plotSpectraChange.md))
  and
  [`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md)/[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md)’s
  `diff_method = "percent_change"` option all promised a percentage but
  never actually multiplied by 100, so values plotted as raw fractions
  (e.g. `-0.057` instead of `-5.7`). Fixed to apply the `*100` when
  percentage output is requested (the default).
- [`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md)/[`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)
  used a simulation’s last *time value* directly as a *row position*,
  which is off by one since time starts at 0 –
  [`plotTotalAbundance()`](https://cmbeese.github.io/mizerReef/reference/plotTotalAbundance.md)
  silently plotted the second-to-last saved timestep, and
  [`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)
  crashed outright for any `mizerReefSim`. Fixed to use the actual last
  row, and switched
  [`plotTotalBiomass()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomass.md)’s
  species subsetting to index by name so it’s robust to the extra
  algae/detritus columns mizerReef adds.
- Fixed several outright crashes and non-functional wrappers, each
  traced to leftover dead code or a copy-paste error:
  [`plotRefugeProfile()`](https://cmbeese.github.io/mizerReef/reference/plotRefugeProfile.md)
  and
  [`plotVulnerable()`](https://cmbeese.github.io/mizerReef/reference/plotVulnerable.md)
  crashed on `MizerSim` input (a later line silently undid an earlier,
  correct branch);
  [`getProductivity()`](https://cmbeese.github.io/mizerReef/reference/getProductivity.md)
  – and therefore
  [`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)/[`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md)
  – never worked on a `MizerSim` object at all
  (`"object 'params' not found"`); and
  [`plotlyProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)/[`plotlyTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md)
  called nonexistent/wrong function names instead of their intended
  [`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)/[`plotTotalBiomassRelative()`](https://cmbeese.github.io/mizerReef/reference/plotTotalBiomassRelative.md)
  counterparts.
- [`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md)’s
  “no species given -\> use all species” default (including the new
  `include_inverts` filtering) used `missing(species)`, which is `FALSE`
  whenever any caller in the chain forwards `species = species` – true
  of every wrapper that calls it – so the “all species” default silently
  never engaged through those wrappers. Switched to `is.null(species)`,
  which reflects intent correctly through any number of wrapper layers.
- [`getProductivity()`](https://cmbeese.github.io/mizerReef/reference/getProductivity.md)‘s
  documented default fishing-length range (7cm minimum, largest species’
  `l_max` maximum) was only applied on the `MizerParams` branch; called
  on a `MizerSim` object (directly, or via
  [`plotProductivity()`](https://cmbeese.github.io/mizerReef/reference/plotProductivity.md),
  [`plot2Productivity()`](https://cmbeese.github.io/mizerReef/reference/plot2Productivity.md),
  [`plotProductivityRelative()`](https://cmbeese.github.io/mizerReef/reference/plotProductivityRelative.md),
  or
  [`plotRelativeContribution()`](https://cmbeese.github.io/mizerReef/reference/plotRelativeContribution.md))
  with `min_fishing_l`/`max_fishing_l` left at their default `NULL`, it
  silently fell through to
  [`mizer::get_size_range_array()`](https://sizespectrum.org/mizer/reference/get_size_range_array.html)’s
  own default of the entire weight spectrum instead, overestimating
  productivity by 30-90% depending on species (confirmed on
  `caribbean_3_model`). Fixed to apply the same default on both
  branches.
- [`plotDegradationScale()`](https://cmbeese.github.io/mizerReef/reference/plotDegradationScale.md)’s
  documentation for a user-supplied custom `trajectory`
  matrix/data.frame described the opposite orientation to what the code
  (and the built-in `rubble_scale`/`algae_scale`/`recovery_scale`
  objects, and `deg_scale` as consumed by
  [`reefDegrade()`](https://cmbeese.github.io/mizerReef/reference/reefDegrade.md))
  actually expect: refuge size bins as rows and bleaching years as
  columns, not bleaching year in the first column. Documentation fixed
  to match actual behaviour; no code change.

#### Compatibility with other mizer extensions (e.g. mizerMR)

- The bundled `caribbean_3_model`/`caribbean_10_model` example models
  now route through mizerReef’s composable rate dispatch when projected
  alongside a second stacked extension, instead of silently falling back
  to a path that assumed no other extension was present. This is what
  caused
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  to error with “incorrect number of dimensions” when combining
  mizerReef with mizerMR.
- [`getSenMort()`](https://cmbeese.github.io/mizerReef/reference/getSenMort.md)
  no longer errors when another extension changes what
  [`initialNResource()`](https://sizespectrum.org/mizer/reference/initialNResource-set.html)
  returns.
- [`getBiomass()`](https://sizespectrum.org/mizer/reference/getBiomass.html)
  for `mizerReefSim` objects now only folds the unstructured algae and
  detritus biomasses into the biomass table, so size-structured
  components added by other extensions (e.g. mizerMR’s multiple-resource
  array) are no longer mis-reshaped into it.
- [`getDegrade()`](https://cmbeese.github.io/mizerReef/reference/getDegrade.md)
  no longer errors when called on a `MizerSim` object.

#### Other

- Fixed
  [`setExtMortParams()`](https://cmbeese.github.io/mizerReef/reference/setExtMortParams.md)
  (and therefore `newReefParams(ext_mort_params = ...)`) so custom
  `ext_mort_params` actually works: it was previously broken for both a
  `data.frame` (crashed on first use) and a named `list` (failed
  validation outright with a misleading error) input.
- **[`reefSenMort()`](https://cmbeese.github.io/mizerReef/reference/reefSenMort.md)/[`getSenMort()`](https://cmbeese.github.io/mizerReef/reference/getSenMort.md)
  computed senescence mortality approximately 5x too high across the
  entire weight spectrum**, for every model, since a since-lost commit
  meant to fix a `NaN`-for-`w`-under-1g edge case. The original (and
  originally documented, thesis-matching) formula is
  `sen_prop * max(0, log10(w)/log10(w_max.i))^sen_curve`, which reaches
  exactly `sen_prop` at a group’s maximum size; the edge-case fix moved
  `sen_prop` inside the floored ratio before it gets raised to the
  `sen_curve` power, so the rate at maximum size became
  `sen_prop^sen_curve` instead (`0.501`, not the intended `0.1`, under
  the package defaults) – a constant `sen_prop^(sen_curve - 1)`
  multiplier at every weight above 1g, not just near maximum size. This
  directly inflated total mortality
  ([`reefMort()`](https://cmbeese.github.io/mizerReef/reference/reefMort.md))
  in every simulation. Fixed to the original formula; the three
  roxygen/vignette copies of this formula (which had independently
  drifted to describe three different, mutually-inconsistent versions,
  none of which matched the buggy code either) are now all consistent
  with it and with each other. `caribbean_3_model` and
  `caribbean_10_model` initially had only
  `other_params$detritus$external` retuned to restore their detritus
  production/consumption balance, which the corrected mortality rate had
  thrown off (senescence mortality contributes to detritus’s “decomp”
  production term) – but that patch alone left fish abundances stale
  relative to the corrected mortality, so neither bundled model was
  actually a genuine steady state of its own current code (confirmed:
  running
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  on them grew biomass 3-9x). Both `caribbean_3_model` and
  `caribbean_10_model` have since been fully recalibrated (fish
  abundances, growth/consumption rates and reproduction parameters all
  re-tuned; verified stable under repeated
  [`reefSteady()`](https://cmbeese.github.io/mizerReef/reference/reefSteady.md)
  calls) – see `inst/scripts/Caribbean_3_model-calibration.R`/
  `inst/scripts/Caribbean_10_model-calibration.R` and
  `vignettes/caribbean_3_model-description.Rmd` for the updated recipes,
  including herbivore `species_params` changes (`satiation`, `age_mat`)
  needed to reach a stable calibration under the corrected mortality.
  `caribbean_10_model` had no surviving calibration script to build from
  (unlike `caribbean_3_model`); the recipe was reconstructed
  interactively, and needed two further fixes beyond
  `caribbean_3_model`’s: (1)
  [`matchReefGrowth()`](https://cmbeese.github.io/mizerReef/reference/matchReefGrowth.md)/[`matchBiomasses()`](https://sizespectrum.org/mizer/reference/matchBiomasses.html)
  had to be applied one species at a time rather than
  mizerExperimental’s usual batched workflow, since batched calls
  destabilised this 10-group model even where they worked fine for
  `caribbean_3_model`’s 3 aggregated groups; and (2) `pred_grab`’s
  stored `age_mat` of 2 years turned out to be inconsistent with its own
  stored von Bertalanffy growth curve and an independently-sourced
  maturity length for its representative species (*Lutjanus apodus*),
  and was corrected to 5.4 years. `eels` and `farm_damsel` – the two
  groups without a FORCE-observed biomass target – were given soft,
  literature-derived targets (not treated as exact) rather than left at
  the near-zero biomass the corrected mortality otherwise drove them to.
  The complexity-increases-biomass relationship central to the
  associated manuscript was explicitly re-verified against both
  recalibrated models, not merely assumed to still hold. See
  `caribbean_10_model`’s documentation and
  `inst/scripts/Caribbean_10_model-calibration.R`’s design note for full
  details and citations.

------------------------------------------------------------------------

See the reference manual and pkgdown site for details on new features
and usage.
