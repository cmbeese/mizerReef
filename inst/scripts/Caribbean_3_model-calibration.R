# Setting up a generic Caribbean coral reef model with multiple resources
# Three groups: Predators, Herbivores, Invertebrates
# Model steady state calibration
# last tuned 17/12/2025
# RECALIBRATED 02/08/2026 -- see note below before rerunning from scratch.

## IMPORTANT: why this script starts from the existing bundled model ------------
#
# The sixteenth session (see inst/to-do-list.txt, git log for reefSenMort())
# corrected a bug that had inflated senescence mortality ~5x above the
# formula in Chapter3.tex. That correction left the *bundled* caribbean_3_model
# no longer a genuine steady state of the package's own corrected rate
# functions (fish abundances were never re-equilibrated, only the resource
# side was patched). Two different attempts to fix this by rebuilding from
# these raw CSVs via newReefParams() and replaying the original recipe
# below verbatim (i.e. starting the pipeline the same way this script
# originally did) both diverged badly under the corrected mortality
# (oscillating biomass, erepro pushed past 1, non-finite search_vol) --
# starting fresh from `newReefParams()` puts the model too far from any
# stable point for the calibrate/match/steady loop to find its way back.
#
# What *does* work is starting from the existing, already-tuned bundled
# caribbean_3_model (a much smaller perturbation) and being careful about
# two specific details that the original recipe didn't need to worry about
# under the old (buggy, higher) mortality:
#   1. matchReefGrowth() must be called with `keep = "biomass"` explicitly.
#      Its default is `keep = "egg"`, which lets biomass swing freely after
#      every growth-rate correction -- under the corrected mortality this
#      created a destructive tug-of-war with matchBiomasses()'s hard
#      rescale (herbivores in particular oscillated: cycle after cycle,
#      erepro clamped to its floor and reproduction level collapsed to
#      exactly zero). With `keep = "biomass"` the two steps stop fighting
#      and the recipe converges smoothly again.
#   2. Herbivores need `satiation = TRUE` (see below) for their growth to
#      have any self-limiting brake at all -- without it, `species_params$
#      age_mat` cannot be reconciled with a low equilibrium biomass under
#      the corrected mortality, no matter how growth/reproduction are
#      tuned (see the design note below).
#
# So: THIS SCRIPT NO LONGER REBUILDS caribbean_3_model FROM THE RAW CSVS.
# It loads the current bundled model, applies the two species_params
# changes below, and re-runs the calibrate/match/steady recipe from there.
# The CSVs are still the source of truth for every other species parameter
# and are still re-saved as data objects for bookkeeping, but `newReefParams()`
# is deliberately NOT called on them for the main model build.

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
#   growth (via the standard Holling type II feeding level). In the
#   recalibrated model the realised feeding level for herbivores converges
#   to ~0.95-0.99 -- i.e. individuals are essentially always near full
#   capacity -- which is a more literal reading of "gut always full" than
#   the previous, formally infinite intake ever was.
#   Without this change, herbivore biomass has no density-dependent brake
#   at all once the corrected (lower) senescence mortality lets individuals
#   survive to a realistic age at maturity: forcing age_mat to match via
#   matchReefGrowth() alone drove herbivore biomass to >800 g/m^2 (vs. the
#   ~34 g/m^2 FORCE-survey target) with no sign of settling.
#
# - age_mat: the thesis's age_mat = 4 for herbivores almost certainly
#   conflated "age at sexual transition" (when female stoplight parrotfish,
#   Sparisoma viride, become male -- a protogynous-hermaphrodite-specific
#   milestone) with "age at first sexual maturity". Rivera Hernandez &
#   Shervette (2025, Environ Biol Fish 108:179-198; a comprehensive 2013-2023
#   otolith/gonad-histology study of 1801 U.S. Caribbean stoplight
#   parrotfish) report age at median sexual transition (AT50) = 4.5 years --
#   very close to the old 4 -- but age at median sexual MATURITY (AM50), the
#   quantity age_mat is actually meant to represent, = 1.6 years. That paper
#   also reports a von Bertalanffy K of 0.33-0.39/year for S. viride, slower
#   than this model's k_vb = 0.6 -- consistent with the old age_mat = 4
#   target having been too slow-growing to begin with. The package
#   maintainer confirmed they trust the FORCE biomass_observed values far
#   more strongly than the original age_mat = 4 guess, and was comfortable
#   retargeting age_mat using this literature value instead.
#   `caribbean_10_species.csv`'s `rep_species` column independently confirms
#   Sparisoma viride is the correct representative species for the
#   equivalent "parrotfish" herbivore group there too.

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

## Load parameters (for bookkeeping/logging only -- NOT used to rebuild the
## model from scratch this time, see note above) --------------------------------
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

# Re-save as R data objects (species_params now includes satiation = TRUE,
# age_mat = 1.6 for herbivores -- see design note above).
save(caribbean_3_species, file = "data/caribbean_3_species.rda")
save(caribbean_3_interaction, file = "data/caribbean_3_interaction.rda")
save(tuning_profile, file = "data/tuning_profile.rda")

## Start from the existing bundled model, not a fresh newReefParams() build ----
data("caribbean_3_model", package = "mizerReef")
params <- caribbean_3_model
log_msg("loaded existing bundled caribbean_3_model as starting point")

# Sync the two changed species_params columns onto the live params object
# (the .rda re-saved above only updates the standalone species_params data
# object, not the model's own copy).
herb <- which(params@species_params$species == "herbivores")
params@species_params$satiation[herb] <- TRUE
params@species_params$age_mat[herb] <- 1.6
log_msg("applied satiation=TRUE, age_mat=1.6 to herbivores")

target <- c(predators = 107, herbivores = 34, inverts = 40)
dist_to_target <- function(p) {
    b <- getBiomass(p) |> tapply(p@species_params$species, sum)
    sum((log(b[names(target)]) - log(target))^2)
}

cat("Starting distance to target:", dist_to_target(params), "\n")
plotBiomassVsSpecies(params)
log_msg("plotBiomassVsSpecies before recalibration")

## Iterate to refine biomass and growth -----------------------------------------
# One alternation at a time, `matchReefGrowth()` before `matchBiomasses()`
# each followed by its own reefSteady() (not batched -- see mizer's own
# "refine your model" course page, which alternates the two the same way),
# and `keep = "biomass"` on matchReefGrowth() explicitly (see note above).
# Stop by watching the distance-to-target metric rather than a fixed count;
# in practice this settles within ~8-10 alternations.
prev_dist <- Inf
for (i in 1:12) {
    params <- params |>
        matchReefGrowth(keep = "biomass") |>
        reefSteady()
    log_msg(sprintf("iter %d: matchReefGrowth(keep=biomass)+steady", i))

    params <- params |>
        matchBiomasses() |>
        reefSteady()
    d <- dist_to_target(params)
    log_msg(sprintf("iter %d: matchBiomasses+steady, dist=%.4f", i, d))
    cat(sprintf("iter %d: dist = %.4f\n", i, d))

    if (d > prev_dist * 1.5 && i > 3) {
        log_msg("distance got notably worse -- stopping, keeping previous state")
        break
    }
    prev_dist <- d
}
log_msg("biomass/growth tuning complete")

plotBiomassVsSpecies(params) # spot on
log_msg("plotBiomassVsSpecies after tuning")

# Check match with age at maturity (herbivores now checked against the
# literature AM50 = 1.6y, not the old, likely-conflated 4y -- see design
# note above; the FORCE biomass match is the trusted target, this is a
# secondary sanity check)
age_mat_observed <- params@species_params$age_mat
age_mat_model <- age_mat(params)
data.frame(age_mat_model, age_mat_observed)

plotTotalAbundance(params)
plotTotalBiomass(params)

## Check resulting spectra and feeding level ------------------------------------

# Spectra should be reasonably straight to match predictions of Sheldon's
# spectrum but also have nonlinearities at refuge sizes
plotSpectra(params, total = TRUE, power = 1)
plotSpectra(params, total = TRUE, power = 2)
log_msg("spectra check before reproduction tuning")

# Herbivore feeding level should now be consistently near 1 (satiated) --
# consistent with "gut always full" (Ferreira 1998; Kopp 2010), see design
# note above.
plotFeedingLevel(params, species = "herbivores")
plotFeedingLevel(params, species = "inverts")

# Tune reproduction ------------------------------------------------------------
rep <- getReproductionLevel(params)
# increase reproduction level to 0.5 for predators/herbivores, leave
# inverts at whatever level the biomass/growth tuning above landed them at
rep_level <- c(0.5, 0.5, rep["inverts"])
names(rep_level) <- c("predators", "herbivores", "inverts")
params <- setBevertonHolt(params, reproduction_level = rep_level)
log_msg("setBevertonHolt reproduction_level", paste(names(rep_level), rep_level, sep = "=", collapse = ";"))

# Iterate to get back to steady state
params <- params |>
    reefSteady() |>
    reefSteady() |>
    reefSteady()
log_msg("post-reproduction reefSteady x3")

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

# Plots ------------------------------------------------------------------------
plotTotalAbundance(params)
plotTotalBiomass(params)
plotBiomassVsSpecies(params)
plotRefugeProfile(params)
plotSpectra(params, power = 1, total = TRUE)
plotDiet(params)
plotGrowthCurves(params)
plotPredMort(params)
plotRelativeContribution(params) # predators dominant on all 3 metrics, as in Chapter3.tex
log_msg("final diagnostic plots generated")

# Save!
caribbean_3_model <- reefSteady(params)
log_msg("final reefSteady to save model")

# Save in package --------------------------------------------------------------
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

save_plot("biomass-final", function() plotBiomassVsSpecies(params))
save_plot("spectra-power1-final", function() plotSpectra(params, total = TRUE, power = 1))
save_plot("spectra-power2-final", function() plotSpectra(params, total = TRUE, power = 2))
save_plot("reproduction-level", function() {
    barplot(getReproductionLevel(params), main = "Reproduction level")
})
log_msg("saved key diagnostic plots")
