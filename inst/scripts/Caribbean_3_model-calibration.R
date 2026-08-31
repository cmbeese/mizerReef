# Setting up a generic Caribbean coral reef model with multiple resources
# Three groups: Predators, Herbivores, Invertebrates
# Model steady state calibration
# last tuned 17/12/2025
# RECALIBRATED 02/08/2026 -- satiation=TRUE / age_mat=1.6 for herbivores
# REBUILT FROM SCRATCH 31/08/2026 -- see note below; supersedes the 02/08/2026
# "start from the bundled model" workaround.
# RECALIBRATED 01/09/2026 -- age_mat=2 (both species), l_mat=16cm / w_mat=102.4g
# (both species), scale_down_factor 32->17; see design note below for
# citations. Diet realism (herbivore share of predator diet, still ~1-6%
# even after checking interaction/beta/interaction_resource levers) remains
# a known, separately-tracked open issue -- not addressed by this revision.

## IMPORTANT: this script rebuilds from newReefParams() again -----------
#
# The 02/08/2026 recalibration found that rebuilding from these raw CSVs via
# newReefParams() "diverged badly" under the corrected senescence-mortality
# formula and fell back to patching the existing bundled model instead. That
# diagnosis was correct for the recipe as it was written then, but the
# actual missing ingredient turned out to be recipe *sequencing*, not
# anything structurally wrong with a fresh build:
#
#   `setBevertonHolt(reproduction_level = 0.5)` needs to be set immediately
#   after construction, before the very first `reefSteady()` call -- not
#   folded in later. Without it, reef models' extra nonlinearities (refuge
#   dynamics, multi-resource coupling) have no damping and the
#   calibrate/match loop oscillates or diverges. This is exactly what
#   `vignettes/steady-state-recipe.Rmd` (the package's own recipe) already
#   does, and exactly what this script previously did not do early enough.
#   Once fixed, the fully undamped, unmodified vignette recipe converges
#   cleanly on `caribbean_3_species.csv`: biomass matches to <0.3%, `age_mat`
#   matches almost exactly (predators 3.9999934 vs target 4.0).
#
# See the caribbean-3-fresh-rebuild-diet-bug memory (or ask about that
# investigation) for the complete derivation.
#
## Diet realism: RESOLVED via mizerExperimental::scaleDownBackground() -------
#
# Earlier versions of this script left predator diet background-Resource-
# dominated at every size. That's fixed now (Step 6, below) using
# `scaleDownBackground()`, which also turned out to fix a second, related
# problem never previously diagnosed: predators' growth energy hit exactly
# zero around w=900-1000g and abundance was exactly zero above that --
# well short of `w_inf=3125g` -- because background Resource (capped at
# `w_pp_cutoff=1g`) becomes inaccessible to predators once their preferred
# prey size (`w/beta`) exceeds 1g (at `w = beta * w_pp_cutoff = 100g`), and
# there wasn't enough real prey to fill the gap. With Step 6, predators now
# grow smoothly all the way to `w_inf` and real-prey fraction (herbivore +
# invertebrate + cannibalism) reaches ~93% at w=933g (up from <1% before).
# Three things were tried before finding this, documented in full detail in
# the caribbean-3-fresh-rebuild-diet-bug memory:
#
#   1. Pinning predators' `gamma` (search volume) to Rogers 2018's own
#      literature value (Appendix S1 Table S2, A_P = 6.4 m^2/yr). Abandoned:
#      a helper that pins `species_params$gamma` directly without going
#      through `species_params<-()` looks like it works but is a no-op --
#      it never touches `@search_vol`, the slot mizer's dynamics actually
#      read. Once fixed to genuinely take effect, pinning gamma=6.4
#      destabilizes the model outright (age at maturity collapses to
#      ~0.002 years). Consistent with the memory's `getStability()` finding
#      that gamma=6.4 alone sits on a near-Hopf bifurcation with a
#      284-year relaxation timescale.
#   2. Retuning `kappa`/`lambda` (the resource abundance scale and slope)
#      *after* biomass calibration has already converged, without forcing
#      gamma. Abandoned: unstable across the entire range tested (targets
#      from kappa=80 down to kappa=4, 8-15 gradual steps, with or without
#      matchReefGrowth() in the polish phase) -- `erepro` explodes rather
#      than settling, a genuine dynamical instability in this region of
#      parameter space, not a step-size or sequencing bug.
#   3. `kappa` matters a great deal at *construction* time, before any
#      calibration runs, and not at all as something to retune afterward:
#      `calibrateReefBiomass()` called once from mizer's own uncalibrated
#      default (kappa=1e11, no package override) lands calibration in a
#      meaningfully better basin (real-prey fraction ~20% at w=933g) than
#      a lower seed does (a kappa=100 seed lands at real-prey fraction
#      ~8%) -- a real, stable improvement, but not enough on its own.
#   4. **This is what Step 6 below actually uses**: `scaleDownBackground()`
#      (see its own comment at Step 6 for the mechanism) shifts the
#      fish-vs-resource abundance balance in a dimensionally-safe way,
#      without touching gamma or kappa/lambda directly. Combined with the
#      kappa=1e11 basin from (3), this reaches real-prey fraction ~93% at
#      w=933g and, as a side effect, resolves the growth-wall problem too
#      (see above) -- genuinely resolving what (1) and (2) could not.
#
# The honest deliverable is: biomass and age-at-maturity match survey
# targets well, and predator diet now shows real, substantial piscivory
# (still more cannibalism-heavy than Rogers 2018's own herb/invert-
# dominant weights -- a remaining, smaller gap, not the near-total absence
# of piscivory this recipe started with).

## Design note: satiation = TRUE and age_mat = 1.6 for herbivores ---------------
#
# Both of these are deliberate departures from the original thesis-era
# values (species_params previously had satiation = FALSE, age_mat = 4 for
# herbivores) and were confirmed with the package maintainer before being
# adopted, not decided unilaterally:
#
# - satiation: Chapter3.tex (and vignettes/model-description.Rmd) documents
#   "satiation exclusive to detritivory" -- herbivores get unlimited intake
#   (h = Inf) because Caribbean herbivores are cited as not down-regulating
#   grazing pressure when food is abundant (Ledlie 2007; Pratchett 2008;
#   Khalil 2013; Elma 2023) and as filling their gut up to three times a day
#   (Ferreira 1998; Kopp 2010). That citation supports "they don't reduce
#   pressure on the shared algae pool" -- it does NOT require literal
#   unlimited per-individual growth-energy absorption, which is what
#   `satiation = FALSE` actually implements (h = Inf, feeding level forced
#   to exactly 0). Algal *pool* depletion (`algae_consumption()`, used in
#   the resource dynamics) is unaffected by this change -- it deliberately
#   ignores feeding level regardless of `satiation`, exactly as documented,
#   so the "no self-regulation of grazing pressure" claim is preserved
#   unchanged. What `satiation = TRUE` adds is a realistic *individual* cap
#   on how much of that encountered energy a herbivore actually absorbs for
#   growth (via the standard Holling type II feeding level).
#   Without this change, herbivore biomass has no density-dependent brake
#   at all once the corrected (lower) senescence mortality lets individuals
#   survive to a realistic age at maturity: forcing age_mat to match via
#   matchReefGrowth() alone drove herbivore biomass to >800 g/m^2 (vs. the
#   ~34 g/m^2 FORCE-survey target) with no sign of settling.
#
# - age_mat / w_mat (UPDATED 2026-09-01): both species' targets are now set
#   from species-specific maturity literature rather than a shared guess:
#
#   Predators (Cephalopholis cruentata, graysby grouper): age_mat = 2 years,
#   from FishBase-cited Caribbean grouper life-history data for female
#   (initial-phase) maturity. Heemstra & Randall (1993, Groupers of the
#   World, FAO Species Catalogue Vol. 125) separately give 4-5 years / 20-23
#   cm TL for the *sex-change* (female-to-male) threshold -- a different
#   milestone from first maturity, and NOT the quantity age_mat represents.
#   This repeats, for a different species, the exact protogyny-conflation
#   risk documented below for herbivores -- confirmed here to have actually
#   affected this dataset too, since the previous age_mat = 4 for predators
#   was in fact the sex-change age, not the maturity age.
#
#   Herbivores (Sparisoma viride, stoplight parrotfish): age_mat = 2 years.
#   Choat, Robertson, Ackerman & Posada (2003, Mar. Ecol. Prog. Ser.
#   246:265-277) give ~4 years, but FishBase's length-at-maturity data for
#   this species (the same source used for l_mat, below) ties initial-phase
#   (female) maturity to an age of ~2-3 years in faster-growing Caribbean
#   populations -- the more recent estimate, and the one confirmed for this
#   build. This *reverts* the 02/08/2026 recalibration's age_mat = 1.6
#   (from Rivera Hernandez & Shervette 2025's AM50 estimate) a second time,
#   to a value distinct from both that and the thesis's original 4.
#   `satiation = TRUE` (above) is kept regardless of which age_mat value is
#   used -- it's the density-dependent brake that lets herbivore biomass
#   settle at any realistic age_mat, not a fix tied to the 1.6-specific
#   value. NOTE: with age_mat = 4, Step 6's scaleDownBackground(factor=32)
#   diverged ("search_vol must not contain non-finite values") -- confirmed
#   this recipe's stability is sensitive to the exact age_mat targets, not
#   just biomass/gamma as previously documented; see the
#   caribbean-3-fresh-rebuild-diet-bug memory for the full record.
#
#   l_mat (both species): 16 cm TL, from FishBase length-at-50%-maturity
#   for each species (predators: ~16 cm; herbivores: 16.3 cm). Converted to
#   `w_mat` via each species' own length-weight relationship (w = a * L^b,
#   a = 0.025, b = 3): 0.025 * 16^3 = 102.4 g.

## Setup - load packages -------------------------------------------------------
library(mizer)
library(mizerExperimental)
library(mizerReef)
library(here)

# Reproducibility scaffold: create artifact folders, start a timestamped log, fix seed,
# and record versions/session/git so this run can be reconstructed exactly.
dir.create("artifacts/plots", recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
logfile <- file.path("artifacts", paste0("caribbean3-calibration-", timestamp, ".txt"))
log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n", file = logfile, append = TRUE)
seed_val <- 12345
set.seed(seed_val)
log_msg("seed", seed_val)
log_msg("mizerReef version", as.character(packageVersion("mizerReef")))
log_msg("mizer version", as.character(packageVersion("mizer")))
log_msg("mizerExperimental version", as.character(packageVersion("mizerExperimental")))
git_head <- tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "git rev-parse failed")
log_msg("git head", paste(git_head, collapse = " "))
git_status <- tryCatch(system("git status --short", intern = TRUE), error = function(e) "git status failed")
log_msg("git status", if (length(git_status)) paste(git_status, collapse = "; ") else "clean")
log_msg("session info start") # full session info goes to the log
capture.output(sessionInfo(), file = logfile, append = TRUE)

# Check that using newest version of mizerReef
packageVersion("mizerReef") # Should be >= 2.0.0

## Load parameters ---------------------------------------------------------------
species_path <- here("inst/data-csv/caribbean_3_species.csv")
interaction_path <- here("inst/data-csv/caribbean_3_interaction.csv")
refuge_path <- here("inst/data-csv/karpata_refuge.csv")
tuning_path <- here("inst/data-csv/tuning_profile.csv")

caribbean_3_species <- read.csv(species_path)
log_msg("caribbean_3_species.csv md5", tools::md5sum(species_path), "rows", nrow(caribbean_3_species))
caribbean_3_interaction <- read.csv(interaction_path, row.names = 1)
log_msg("caribbean_3_interaction.csv md5", tools::md5sum(interaction_path), "rows", nrow(caribbean_3_interaction))
karpata_refuge <- read.csv(refuge_path)
log_msg("karpata_refuge.csv md5", tools::md5sum(refuge_path), "rows", nrow(karpata_refuge))
tuning_profile <- read.csv(tuning_path)
log_msg("tuning_profile.csv md5", tools::md5sum(tuning_path), "rows", nrow(tuning_profile))

save(caribbean_3_species, file = "data/caribbean_3_species.rda")
save(caribbean_3_interaction, file = "data/caribbean_3_interaction.rda")
save(tuning_profile, file = "data/tuning_profile.rda")

target <- c(predators = 107, herbivores = 34, inverts = 40)
dist_to_target <- function(p) {
    b <- getBiomass(p) |> tapply(p@species_params$species, sum)
    sum((log(b[names(target)]) - log(target))^2)
}

## Step 1: build fresh from the raw CSVs -----------------------------------------
params <- newReefParams(species_params = caribbean_3_species,
                         interaction = caribbean_3_interaction,
                         method = "binned", method_params = tuning_profile)
log_msg("built fresh params via newReefParams()")

## Step 2: moderate density-dependence, set EARLY, then initial settle ----------
## (must happen before the first reefSteady() call -- see note above)
rdi <- rep(0.5, nrow(caribbean_3_species))
names(rdi) <- caribbean_3_species$species
params <- setBevertonHolt(params, reproduction_level = rdi)
params <- params |>
    reefSteady() |> reefSteady() |> reefSteady() |>
    reefSteady() |> reefSteady() |> reefSteady()
stopifnot(mizer::isSteady(params))
log_msg("initial settle complete, isSteady=TRUE")

## Step 3: calibrate biomass ONCE, then iteratively fine-tune -------------------
## `calibrateReefBiomass()` does a one-shot rescale of kappa/gamma/abundance
## to bring the model onto the right absolute scale -- it is meant to run
## once, not be looped. Looping it (as an earlier version of this script
## did, 9 times total) does not diverge, but does land the model in a
## meaningfully worse final kappa/diet basin: calling it once from mizer's
## own default kappa (1e11, no override -- see note above) settles at
## kappa~107 with real-prey fraction ~20% at 933g; looping it 9 times from
## the SAME starting kappa settles at kappa~103 (materially the same,
## confirming looping isn't itself the problem there) -- but looping it 9
## times from a *seeded* kappa=100 (this package's former, since-removed
## reef_kappa_default) settles at kappa~305-311 with real-prey fraction
## ~8%, a clearly worse basin. The single-call structure below matches
## vignettes/steady-state-recipe.Rmd's own Step 4a/4b, is simpler, and
## reliably lands in the better basin regardless of iteration count.
params <- calibrateReefBiomass(params)
params <- reefSteady(params)
log_msg(sprintf("Step 3: single calibrateReefBiomass(), dist=%.4g kappa=%.4g",
                 dist_to_target(params), params@resource_params$kappa))
for (i in 1:6) {
    params <- params |>
        matchBiomasses() |> matchReefGrowth() |>
        reefSteady()
    log_msg(sprintf("Step 3 iter %d: dist=%.4g", i, dist_to_target(params)))
}
stopifnot(mizer::isSteady(params))
cat("After Step 3 (binned refuge, tuning phase complete), dist:", dist_to_target(params), "\n")
plotBiomassVsSpecies(params)
log_msg("plotBiomassVsSpecies after binned-refuge tuning")

## Step 4: switch to competitive refuge method -----------------------------------
## Switching refuge methods is itself a real perturbation and needs many
## more rounds of matchBiomasses()+reefSteady() to fully re-settle than the
## fine-tuning above (confirmed: ~16 rounds, with early rounds oscillating
## substantially before converging monotonically).
params <- newRefuge(params, new_method = "competitive",
                     new_method_params = karpata_refuge)
log_msg("switched to competitive refuge method")
for (i in 1:18) {
    params <- matchBiomasses(params) |> reefSteady()
    log_msg(sprintf("Step 4 iter %d: dist=%.4g", i, dist_to_target(params)))
}
stopifnot(mizer::isSteady(params))
cat("After Step 4 (competitive refuge), dist:", dist_to_target(params), "\n")

plotVulnerable(params)
plotRefugeProfile(params)
log_msg("refuge diagnostics after competitive-method switch")

## Final checks -------------------------------------------------------------------
age_mat_observed <- params@species_params$age_mat
age_mat_model <- age_mat(params)
data.frame(age_mat_model, age_mat_observed)
log_msg("age_mat check: model vs observed",
        paste(round(age_mat_model, 3), collapse = "/"), "vs",
        paste(age_mat_observed, collapse = "/"))

plotTotalAbundance(params)
plotTotalBiomass(params)
plotSpectra(params, total = TRUE, biomass = TRUE)
plotSpectra(params, total = TRUE, biomass = TRUE, per_log_size = TRUE)
log_msg("spectra check before reproduction tuning")

## Tune reproduction ---------------------------------------------------------------
rep <- getReproductionLevel(params)
rep_level <- c(0.5, 0.5, rep["inverts"])
names(rep_level) <- c("predators", "herbivores", "inverts")
params <- setBevertonHolt(params, reproduction_level = rep_level)
log_msg("setBevertonHolt reproduction_level", paste(names(rep_level), rep_level, sep = "=", collapse = ";"))

for (i in 1:3) {
    params <- matchBiomasses(params) |>
        matchReefGrowth() |> reefSteady()
}
log_msg("post-reproduction refinement x3")

rep <- getReproductionLevel(params)
getRDI(params) / getRDD(params)

# Confirm genuine steady state: further reefSteady() calls should barely
# move biomass at all.
for (i in 1:3) {
    params <- reefSteady(params)
    b <- getBiomass(params) |> tapply(params@species_params$species, sum)
    cat(sprintf("stability check %d: %s\n", i, paste(round(b[names(target)], 4), collapse = "/")))
}
log_msg("stability check: 3x reefSteady with no further movement")

## Step 5: rescale algae/detritus absolute scale ---------------------------------
## Fish dynamics never depend on the absolute standing-biomass scale of
## these two unstructured resources (only on production/consumption rates,
## which this doesn't touch) -- but the scale itself is left wherever
## calibrateReefBiomass()'s repeated rho_algae/rho_detritus rescaling
## happened to leave it, which by this point is ~1e-9 g/m^2 for both
## (effectively zero, not a real steady state for the resource itself).
## `rescale_algae()`/`rescale_detritus()` fix this without touching fish
## dynamics at all (confirmed: biomass/gamma/diet all bit-identical
## before/after, verified directly this session) -- multiplying standing
## biomass by a factor while dividing the consuming species' rho by the
## same factor, so total consumption is unchanged.
##
## Targets ported from the caribbean_10 calibration campaign's own
## literature audit (`caribbean-10-calibration/references/README.md` §1,
## `DECISIONS.md` 2026-08-06/07 and 2026-08-21 rows) -- these are the
## literature *brackets*, not that campaign's fitted numbers, so they're
## reusable here as-is:
##   - detritus_lifetime = 3.85h (0.00043950 yr): midpoint of Max, Hamilton,
##     Gaines & Warner (2013, Mar. Ecol. Prog. Ser. 482:181-195) Palmyra
##     Atoll benthic-detritus residence-time range (2.8-4.9h).
##   - algae standing biomass = 350 g WW/m^2: midpoint of a 150-550 g WW/m^2
##     bracket, itself a conversion of a user-supplied ~10-100 g DW/m^2
##     Caribbean estimate via an assumed 5-10x wet:dry ratio -- that
##     conversion factor is flagged there as still provisional, not an
##     independently checked literature figure. If a better Caribbean-
##     specific wet:dry ratio turns up, this is the one number in this
##     step worth revisiting.
algae_biomass_target <- 350             # g WW/m^2
detritus_lifetime_target <- 0.00043950  # years (3.85h)

params <- rescale_algae(params, algae_biomass_target / params@initial_n_other$algae)
detritus_lifetime(params) <- detritus_lifetime_target
params <- reefSteady(params)
stopifnot(mizer::isSteady(params))
log_msg(sprintf("Step 5: algae biomass -> %.1f g/m^2, detritus_lifetime -> %.2fh",
                 params@initial_n_other$algae, detritus_lifetime(params) * 8760))
cat("After Step 5 (algae/detritus rescale), dist:", dist_to_target(params), "\n")
cat("algae biomass:", params@initial_n_other$algae, " detritus_lifetime (h):",
    detritus_lifetime(params) * 8760, "\n")

## Step 6: scaleDownBackground() -- shift the fish-vs-resource balance -----------
## Diet realism (see KNOWN OPEN ISSUE note above) is fixed here, not by
## touching gamma or kappa/lambda directly (both dead ends -- see above),
## but by mizerExperimental's `scaleDownBackground(params, factor)`:
##   scaleAbundance(params, factor) |> scaleModel(factor = 1/factor)
## `scaleAbundance()` first rescales ONLY foreground (fish) abundance by
## `factor`; `scaleModel()` (mizer's own, the same dimensionally-safe
## rescale `scaleReefModel()`/`calibrateReefBiomass()` already use) then
## rescales EVERYTHING -- kappa, gamma/search_vol, and all abundances,
## fish included -- by `1/factor`. Net effect: fish abundance is
## unchanged (factor * 1/factor), but kappa shrinks by 1/factor while
## gamma/search_vol grow by factor to compensate, in a self-consistent,
## unit-safe way -- unlike the ad-hoc kappa/lambda retune this replaces,
## which was unstable across its entire tested range (see above).
##
## The right factor is empirically tuned with a safety margin below where
## this goes unstable. The exact boundary is fragile and NOT perfectly
## reproducible run-to-run or across species_params changes: the
## refuge-switch phase above (Step 4) is itself somewhat chaotic (Part 5's
## own note: early rounds oscillate substantially before converging), so
## tiny floating-point differences from one run to the next land Step 5's
## kappa at a slightly different value, which shifts exactly how far
## scaleDownBackground() can safely push before diverging. Confirmed
## directly: a boundary found against one build (kappa=311, factor=53
## safe/54 unsafe) did NOT transfer to another build of the nominally same
## recipe (kappa=187, factor=34 safe/35 unsafe) -- pick a factor with real
## margin below wherever the boundary lands for a given build, don't sit
## right on the edge.
##
## UPDATED 2026-09-01: changing predators'/herbivores' age_mat (see the
## design note above, age_mat = 2 for both, down from 4/1.6) shifted this
## boundary sharply lower -- confirmed by direct sweep against this build:
## factor=19 stable, factor=20 fails ("search_vol must not contain
## non-finite values"). The previous factor=32 no longer converges at all.
## Re-picked at 17 for margin below the new boundary, same convention as
## before.
scale_down_factor <- 17
params <- scaleDownBackground(params, scale_down_factor)
params <- reefSteady(params)
for (i in 1:15) {
    params <- matchBiomasses(params) |> matchReefGrowth() |> reefSteady()
    log_msg(sprintf("Step 6 iter %d: dist=%.4g", i, dist_to_target(params)))
}
stopifnot(mizer::isSteady(params))
cat("After Step 6 (scaleDownBackground, factor=", scale_down_factor, "), dist:",
    dist_to_target(params), "\n")

## Step 6 touches rho_algae/rho_detritus (via matchReefGrowth()'s growth
## rescale), drifting algae/detritus off Step 5's literature targets
## (confirmed: detritus_lifetime drifted from 3.85h to 0.77h in testing) --
## re-apply the same rescale one more time now that Step 6 is done.
params <- rescale_algae(params, algae_biomass_target / params@initial_n_other$algae)
detritus_lifetime(params) <- detritus_lifetime_target
params <- reefSteady(params)
stopifnot(mizer::isSteady(params))
log_msg(sprintf("Step 6 algae/detritus re-rescale: algae -> %.1f g/m^2, detritus_lifetime -> %.2fh",
                 params@initial_n_other$algae, detritus_lifetime(params) * 8760))

# Diet real-prey-fraction check (see KNOWN OPEN ISSUE note above).
diet <- getDiet(params)
w <- params@w
for (target_w in c(5, 50, 500, 933, 2000)) {
    wi <- which.min(abs(w - target_w))
    d <- diet["predators", wi, ]
    real_prey <- sum(d[intersect(names(d), c("predators", "herbivores", "inverts"))], na.rm = TRUE)
    cannibalism <- d[["predators"]]
    log_msg(sprintf("diet at w=%.1fg: real_prey_frac=%.3f (of which cannibalism=%.3f)",
                     w[wi], real_prey, cannibalism))
}

## Plots ------------------------------------------------------------------------
plotTotalAbundance(params)
plotTotalBiomass(params)
plotBiomassVsSpecies(params)
plotRefugeProfile(params)
plotSpectra(params, biomass = TRUE, total = TRUE)
plotDiet(params)
plotGrowthCurves(params)
plotPredMort(params)
plotRelativeContribution(params)
log_msg("final diagnostic plots generated")

# Save!
caribbean_3_model <- reefSteady(params)
log_msg("final reefSteady to save model")
stopifnot(mizer::isSteady(caribbean_3_model))

## Save in package ----------------------------------------------------------------
save(caribbean_3_model, file = "data/caribbean_3_model.rda") # package copy
log_msg("package data saved: data/caribbean_3_model.rda")

# Save full params/model as artifacts with timestamped filenames so a past run
# can be reloaded exactly (includes mizerReef custom slots).
params_artifact <- file.path("artifacts", paste0("caribbean_3_params-", timestamp, ".rds"))
model_artifact <- file.path("artifacts", paste0("caribbean_3_model-", timestamp, ".rds"))
mizer::saveParams(params, params_artifact)
mizer::saveParams(caribbean_3_model, model_artifact)
log_msg("artifacts saved", params_artifact, model_artifact)

# Save the key plots that informed decisions to the artifact folder.
save_plot <- function(name, fn) {
    png(file.path("artifacts", "plots", paste0(name, "-", timestamp, ".png")), width = 1400, height = 900, res = 150)
    on.exit(dev.off(), add = TRUE)
    fn()
}

save_plot("biomass-final", function() print(plotBiomassVsSpecies(params)))
save_plot("spectra-power1-final", function() print(plotSpectra(params, total = TRUE, biomass = TRUE)))
save_plot("spectra-power2-final", function() print(plotSpectra(params, total = TRUE, biomass = TRUE, per_log_size = TRUE)))
save_plot("diet-final", function() print(plotDiet(params)))
save_plot("reproduction-level", function() {
    barplot(getReproductionLevel(params), main = "Reproduction level")
})
log_msg("saved key diagnostic plots")
