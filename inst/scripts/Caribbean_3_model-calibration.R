# Setting up a generic Caribbean coral reef model with multiple resources
# Three groups: Predators, Herbivores, Invertebrates
# Model steady state calibration
# last tuned 17/12/2025

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

## Load parameters -------------------------------------------------------------

# Load species parameter data and save as an object to be included with package.
# We log md5 checksums and row counts to prove which input files were used.
species_path <- here("inst/data-csv/caribbean_3_species.csv")
interaction_path <- here("inst/data-csv/caribbean_3_interaction.csv")
refuge_path <- here("inst/data-csv/karpata_refuge.csv")
tuning_path <- here("inst/data-csv/tuning_profile.csv")

caribbean_3_species <- read.csv(species_path)
log_msg("caribbean_3_species.csv md5", tools::md5sum(species_path), "rows", nrow(caribbean_3_species))
caribbean_3_interaction <- read.csv(interaction_path, row.names = 1)
log_msg("caribbean_3_interaction.csv md5", tools::md5sum(interaction_path), "rows", nrow(caribbean_3_interaction))
# Refuge densities from Karpata reef in Bonaire, FORCE dataset
karpata_refuge <- read.csv(refuge_path)
log_msg("karpata_refuge.csv md5", tools::md5sum(refuge_path), "rows", nrow(karpata_refuge))
tuning_profile <- read.csv(tuning_path) # 60% refuge for all size classes
log_msg("tuning_profile.csv md5", tools::md5sum(tuning_path), "rows", nrow(tuning_profile))

# Save these as R data objects
save(caribbean_3_species, file = "data/caribbean_3_species.rda")
save(caribbean_3_interaction, file = "data/caribbean_3_interaction.rda")
save(karpata_refuge, file = "data/caribbean_3_refuge.rda")
save(tuning_profile, file = "data/tuning_profile.rda")

# With these parameters, herbivores consume plankton at small sizes and
#   transition fully to algae by maturity
# With these parameters, invertebrates consume plankton and detritus,
#   with the proportion of detritus increasing with size

## Set model -------------------------------------------------------------------
params <- newReefParams(
    species_params = caribbean_3_species,
    interaction = caribbean_3_interaction,
    method = "binned",
    method_params = tuning_profile
)
log_msg("newReefParams set with binned method")

## Project to first steady state -----------------------------------------------
params <- reefSteady(params)
log_msg("reefSteady initial complete")

## Calibrate biomasses and growth ----------------------------------------------

# Match observed species group biomasses
params <- calibrateReefBiomass(params)
params <- matchBiomasses(params)
params <- reefSteady(params)
log_msg("calibrateReefBiomass + matchBiomasses cycle 1 complete")

# Match observed growth rates
params <- matchReefGrowth(params)
params <- reefSteady(params)
log_msg("matchReefGrowth cycle 1 complete")

# Check for match with age at maturity
age_mat_observed <- caribbean_3_species$age_mat
age_mat_model <- age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Not bad

# Check biomass match
plotBiomassVsSpecies(params)
log_msg("plotBiomassVsSpecies initial check")
# Biomasses way off from observations

# Iterate to refine biomass - repeat many times; piping keeps it readable.
params <- params |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady()
log_msg("biomass/growth tuning (binned method) complete")

plotBiomassVsSpecies(params) # spot on
log_msg("plotBiomassVsSpecies after binned tuning")

# Check match with observed age at maturity
age_mat_observed <- caribbean_3_species$age_mat
age_mat_model <- age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Closer than needed

plotTotalAbundance(params)
plotTotalBiomass(params)

## Now switch to competitive method --------------------------------------------
params <- newRefuge(params,
    new_method = "competitive",
    new_method_params = caribbean_3_refuge
)
log_msg("switched refuge to competitive")

# Match biomasses again - run three times
params <- params |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady() |>
    calibrateReefBiomass() |>
    matchBiomasses() |>
    matchReefGrowth() |>
    reefSteady()
log_msg("biomass/growth tuning (competitive refuge) complete")

# Make sure new refuge is in place
plotVulnerable(params)
plotRefugeProfile(params)
log_msg("checked refuge/vulnerable profiles")
# looks good

plotBiomassVsSpecies(params) # spot on
log_msg("plotBiomassVsSpecies after competitive tuning")

# Check match with observed age at maturity
age_mat_observed <- caribbean_3_species$age_mat
age_mat_model <- age_mat(params)
data.frame(age_mat_model, age_mat_observed)
# Also spot on

## Check resulting spectra and tune resources ----------------------------------

# Spectra should be reasonably straight to match predictions of Sheldon's
# spectrum but also have nonlinearities at refuge sizes
plotSpectra(params, total = TRUE, power = 1) # looks straight, some bumps
plotSpectra(params, total = TRUE, power = 2) # resource looks low
log_msg("spectra check before reproduction tuning")

# plot feeding level to check if resource is too low
plotFeedingLevel(params, species = "inverts")
# Invertebrate feeding level is stable throughout life - there is enough
# resource to support fish, not too little or too much

# Tune reproduction ------------------------------------------------------------
# We do not have yield or catch data - can't tune size distribution
# First attempt to set very low to see what the minimum values are
params <- setBevertonHolt(params, erepro = 0.0001)
# Now set setting erepro same for all species, as low as possible
params <- setBevertonHolt(params, erepro = 0.0134)
log_msg("setBevertonHolt erepro sequence 0.0001 -> 0.0134")
params <- reefSteady(params)
# Check reproduction level (value between 0 and 1) - should be higher for
# larger, slow growing species and low for small, fast growing ones
rep <- getReproductionLevel(params)
# These are very low for predators, should be higher and too high
# for invertebrates. A reproduction level closer to one means reproduction
# rate is almost totally independent of the investment into reproduction
# Reproduction should be somewhat density independent on reefs

# Check comparison of density dependent & independent reduction
getRDI(params) / getRDD(params)
# preds, 1:1 - maybe too density dependent
# inverts, herbs 20+:1 more density independent - maybe too much

# increase reproduction level to 0.5 for all
rep_level <- c(0.5, 0.5, rep["inverts"])
names(rep_level) <- c("predators", "herbivores", "inverts")
params <- setBevertonHolt(params,
    reproduction_level = rep_level
)
log_msg("setBevertonHolt reproduction_level", paste(names(rep_level), rep_level, sep = "=", collapse = ";"))

# Iterate to get back to steady state
params <- params |>
    reefSteady() |>
    reefSteady() |>
    reefSteady()
log_msg("post-reproduction reefSteady x3")

# Check new reproduction - these look better
rep <- getReproductionLevel(params)
getRDI(params) / getRDD(params)
# Now density independent reproduction is double for predators
# And reproduction is somewhat density dependent for herbs and inverts

# Check new spectra
plotSpectra(params, total = TRUE, power = 1)
plotSpectra(params, total = TRUE, power = 2)
log_msg("spectra check after reproduction tuning")
# These still look good

# Plots ------------------------------------------------------------------------
plotTotalAbundance(params) # Total abundances look reasonable, inverts in range
plotTotalBiomass(params)
plotBiomassVsSpecies(params)
plotRefuge(params)
plotSpectra(params, power = 1, total = TRUE)
plotDiet(params)
plotGrowthCurves(params)
plotPredMort(params)
log_msg("final diagnostic plots generated")
# Everything looks good here! I am happy with my results.

# Save!
caribbean_3_model <- reefSteady(params)
log_msg("final reefSteady to save model")

# Save in package --------------------------------------------------------------
#
# caribbean_3_model <- caribbean_3_model
# caribbean_3_model@other_params$degrade <- FALSE
# caribbean_3_model@other_params$new_refuge <- FALSE

# Params object
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
