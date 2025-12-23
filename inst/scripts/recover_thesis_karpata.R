# Recover legacy thesis-karpata model saved without custom slots
# 
# PURPOSE:
# This model was calibrated with an older version of mizerReef (from main branch)
# before custom slots were added (major_update branch). Loading it directly with
# the current package fails due to class definition mismatches. This script bridges
# the gap by loading it with the matching old package, extracting all parameters
# and diagnostics, then reconstructing it under the current package.
#
# WORKFLOW:
# 1) Install legacy mizerReef (main branch) into a temporary library
# 2) Load the old .rda file with the matching class definition
# 3) Extract and document ALL parameters, slots, diagnostics, and plots
# 4) Save portable data (species params, interaction, etc.) to neutral RDS
# 5) Reload under current mizerReef (major_update) and rebuild the model
# 6) Generate comparison plots to identify calibration differences
#
# OUTPUT:
# - artifacts/legacy-diagnostics-thesis-karpata.txt: Full parameter dump
# - artifacts/legacy-plots/: All diagnostic plots from the legacy model
# - artifacts/legacy_dump-thesis-karpata.rds: Raw portable data
# - artifacts/thesis-karpata-modern.rds/.rda: Modernized version
#
# USE CASE:
# When you can't remember which calibration "test" was saved, inspect the
# diagnostics file and plots to identify parameter choices (erepro, resource
# scaling, refuge method, etc.) and use them as starting points for new runs.

suppressPackageStartupMessages({
    library(withr)   # for with_libpaths to isolate temp library
    library(remotes) # for install_github
    library(here)    # for portable paths
})

# Null-coalescing operator: returns y if x is NULL, otherwise x
`%||%` <- function(x, y) if (is.null(x)) y else x

# Define file paths for all inputs and outputs
legacy_rda <- here("inst/archive/data-archive/thesis-karpata_model.rda")  # Source: old model
legacy_dump_rds <- here("artifacts/legacy_dump-thesis-karpata.rds")       # Portable data dump
legacy_diagnostics_txt <- here("artifacts/legacy-diagnostics-thesis-karpata.txt")  # Full report
modern_params_rds <- here("artifacts/thesis-karpata-modern.rds")          # Modernized .rds
modern_params_rda <- here("artifacts/thesis-karpata-modern.rda")          # Modernized .rda

# Create output directories if they don't exist
dir.create(here("artifacts"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("artifacts/legacy-plots"), recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# STEP 1: LOAD LEGACY MODEL WITH MATCHING PACKAGE VERSION
# ============================================================================
# We need to install the old mizerReef (main branch) into a temporary library
# so we can load the old .rda file without class mismatch errors. We isolate
# this in a temp lib to avoid contaminating the current R session.

legacy_lib <- tempfile("legacy_lib_")  # Create temp directory for legacy package
dir.create(legacy_lib)

cat("=== LEGACY MODEL RECOVERY ===\n")
cat("Installing legacy mizerReef from main branch...\n")

# with_libpaths temporarily adds legacy_lib to the library search path
# so packages installed within this block go there and don't affect your main library
with_libpaths(new = legacy_lib, action = "prefix", code = {
    
    # Install mizerReef from main branch (the version before custom slots were added)
    remotes::install_github("cmbeese/mizerReef@main", upgrade = "never", quiet = FALSE)
    
    # Load the legacy packages now that they're installed in the temp lib
    library(mizer)
    library(mizerReef)
    
    cat("\nLoading legacy model from:", legacy_rda, "\n")
    obj_name <- load(legacy_rda)  # Returns the name of the loaded object
    legacy <- get(obj_name)        # Get the actual object by name
    
    cat("Successfully loaded object:", obj_name, "\n")
    cat("Class:", class(legacy), "\n\n")
    
    # ========================================================================
    # EXTRACT EVERYTHING FOR FORENSIC ANALYSIS
    # ========================================================================
    # We need to capture all parameter values, slots, and diagnostics while
    # we still have the object loaded with the matching class definition.
    # Once we exit this with_libpaths block, we won't be able to load it again.
    
    sink(legacy_diagnostics_txt)  # Redirect all output to text file
    cat("=== LEGACY MODEL DIAGNOSTIC REPORT ===\n")
    cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
    cat("Source:", legacy_rda, "\n\n")
    
    # Record which package versions were used to create this model
    cat("--- PACKAGE VERSIONS (legacy environment) ---\n")
    cat("mizerReef:", as.character(packageVersion("mizerReef")), "\n")
    cat("mizer:", as.character(packageVersion("mizer")), "\n\n")
    
    # List all slots in the object (reveals structure differences)
    cat("--- SLOT NAMES ---\n")
    print(slotNames(legacy))
    cat("\n")
    
    # Dump all species parameters (includes erepro, reproduction_level if stored)
    cat("--- SPECIES PARAMETERS ---\n")
    print(legacy@species_params)
    cat("\n")
    
    # Interaction matrix (who eats whom)
    cat("--- INTERACTION MATRIX ---\n")
    print(legacy@interaction)
    cat("\n")
    
    # Resource parameters (background spectrum, carrying capacity)
    cat("--- RESOURCE PARAMETERS ---\n")
    print(legacy@resource_params)
    cat("\n")
    
    # other_params contains custom settings: refuge method, degrade flags, etc.
    cat("--- OTHER PARAMETERS ---\n")
    print(str(legacy@other_params))
    cat("\n")
    
    # ========================================================================
    # EXTRACT KEY DIAGNOSTICS
    # ========================================================================
    # Try to extract reproduction, mortality, biomass diagnostics.
    # Wrapped in tryCatch in case functions don't exist in legacy version.
    
    tryCatch({
        cat("--- REPRODUCTION LEVEL ---\n")
        print(getReproductionLevel(legacy))
        cat("\n")
    }, error = function(e) cat("Could not extract reproduction level:", e$message, "\n\n"))
    
    tryCatch({
        cat("--- RDI / RDD RATIO ---\n")
        print(getRDI(legacy) / getRDD(legacy))
        cat("\n")
    }, error = function(e) cat("Could not extract RDI/RDD:", e$message, "\n\n"))
    
    tryCatch({
        cat("--- AGE AT MATURITY ---\n")
        print(age_mat(legacy))
        cat("\n")
    }, error = function(e) cat("Could not extract age at maturity:", e$message, "\n\n"))
    
    cat("--- SPECIES BIOMASSES ---\n")
    tryCatch({
        print(getBiomass(legacy))
        cat("\n")
    }, error = function(e) cat("Could not extract biomass:", e$message, "\n\n"))
    
    sink()  # Stop redirecting output
    cat("Diagnostics written to:", legacy_diagnostics_txt, "\n")
    
    # ========================================================================
    # SAVE DIAGNOSTIC PLOTS
    # ========================================================================
    # Generate and save all key plots that you would have inspected during
    # calibration. Compare these to your memory of "test 2" to identify which
    # test the legacy model represents.
    
    cat("Generating legacy plots...\n")
    
    # Biomass by species (should match observed values if calibration worked)
    tryCatch({
        png(here("artifacts/legacy-plots/biomass-vs-species.png"), width = 1400, height = 900, res = 150)
        print(plotBiomassVsSpecies(legacy))
        dev.off()
    }, error = function(e) cat("Could not plot biomass:", e$message, "\n"))
    
    # Spectra plots (should follow Sheldon's spectrum, reveal resource levels)
    tryCatch({
        png(here("artifacts/legacy-plots/spectra-power1.png"), width = 1400, height = 900, res = 150)
        print(plotSpectra(legacy, total = TRUE, power = 1))
        dev.off()
    }, error = function(e) cat("Could not plot spectra power 1:", e$message, "\n"))
    
    tryCatch({
        png(here("artifacts/legacy-plots/spectra-power2.png"), width = 1400, height = 900, res = 150)
        print(plotSpectra(legacy, total = TRUE, power = 2))
        dev.off()
    }, error = function(e) cat("Could not plot spectra power 2:", e$message, "\n"))
    
    # Feeding level (indicates if resources are sufficient)
    tryCatch({
        png(here("artifacts/legacy-plots/feeding-level.png"), width = 1400, height = 900, res = 150)
        print(plotFeedingLevel(legacy))
        dev.off()
    }, error = function(e) cat("Could not plot feeding level:", e$message, "\n"))
    
    # Diet composition (what each species eats)
    tryCatch({
        png(here("artifacts/legacy-plots/diet.png"), width = 1400, height = 900, res = 150)
        print(plotDiet(legacy))
        dev.off()
    }, error = function(e) cat("Could not plot diet:", e$message, "\n"))
    
    # Predation mortality (death rates by size)
    tryCatch({
        png(here("artifacts/legacy-plots/pred-mort.png"), width = 1400, height = 900, res = 150)
        print(plotPredMort(legacy))
        dev.off()
    }, error = function(e) cat("Could not plot pred mort:", e$message, "\n"))
    
    # Vulnerable abundance (what's exposed to predation after refuge)
    tryCatch({
        png(here("artifacts/legacy-plots/vulnerable.png"), width = 1400, height = 900, res = 150)
        print(plotVulnerable(legacy))
        dev.off()
    }, error = function(e) cat("Could not plot vulnerable:", e$message, "\n"))
    
    # Refuge profile (density by size class in refuge)
    tryCatch({
        png(here("artifacts/legacy-plots/refuge-profile.png"), width = 1400, height = 900, res = 150)
        print(plotRefugeProfile(legacy))
        dev.off()
    }, error = function(e) cat("Could not plot refuge profile:", e$message, "\n"))
    
    cat("Legacy plots saved to artifacts/legacy-plots/\n\n")
    
    # ========================================================================
    # DUMP PORTABLE DATA
    # ========================================================================
    # Extract all the core data structures that don't depend on the class
    # definition. We can reload these under the current package without errors.
    # This includes all the parameters that define the model state.
    
    legacy_dump <- list(
        species_params  = legacy@species_params,      # Species traits, erepro, etc.
        interaction     = legacy@interaction,          # Predation interaction matrix
        gear_params     = legacy@gear_params,          # Fishing gear parameters
        resource_params = legacy@resource_params,      # Background resource spectrum
        other_params    = legacy@other_params,         # Custom mizerReef settings
        w_min_idx       = legacy@w_min_idx,            # Minimum size indices per species
        w_pp_cutoff_idx = legacy@w_pp_cutoff_idx,      # Resource cutoff index
        # Include raw slot data that might have clues about package version
        metadata = list(
            slot_names = slotNames(legacy),  # All slot names in legacy object
            class = class(legacy)             # Class hierarchy
        )
    )
    saveRDS(legacy_dump, legacy_dump_rds)
    cat("Legacy dump written to:", legacy_dump_rds, "\n\n")
})


# ============================================================================
# STEP 2: REHYDRATE UNDER CURRENT MIZERREEF (MAJOR_UPDATE BRANCH)
# ============================================================================
# Now that we've extracted all the legacy data, we reload it and rebuild
# the model using the current package version with the new custom slots.
# This creates a modern version that can be used with current code.

cat("=== REHYDRATING UNDER CURRENT MIZERREEF ===\n")
suppressPackageStartupMessages({
    library(mizer)
    library(mizerExperimental)
    library(mizerReef)
})

cat("Current package versions:\n")
cat("mizerReef:", as.character(packageVersion("mizerReef")), "\n")
cat("mizer:", as.character(packageVersion("mizer")), "\n\n")

# Load the portable dump we created in step 1
legacy <- readRDS(legacy_dump_rds)

cat("Reconstructing params from legacy dump...\n")

# ========================================================================
# INFER CALIBRATION SETTINGS FROM LEGACY DATA
# ========================================================================
# Try to figure out which refuge method and parameters were used.
# These are stored in other_params and may have different names in legacy.

method_guess <- legacy$other_params$new_refuge_method %||% "binned"
method_params_guess <- legacy$other_params$tuning_profile %||%
    legacy$other_params$refuge_profile %||%
    legacy$other_params$method_params %||% NULL

cat("Detected refuge method:", method_guess, "\n")

# ========================================================================
# REBUILD MODEL WITH CURRENT PACKAGE
# ========================================================================
# newReefParams creates a new MizerReefParams object using the current
# class definition, populated with the legacy parameter values.

params_new <- newReefParams(
    species_params = legacy$species_params,
    interaction = legacy$interaction,
    method = method_guess,
    method_params = method_params_guess
)

# ========================================================================
# RESTORE LEGACY FLAGS/SETTINGS
# ========================================================================
# Copy over any custom flags that were set in the legacy version.
# These control special behaviors like degradation or refuge handling.

if (!is.null(legacy$other_params$degrade)) {
    params_new@other_params$degrade <- legacy$other_params$degrade
    cat("Set degrade flag to:", legacy$other_params$degrade, "\n")
}
if (!is.null(legacy$other_params$new_refuge)) {
    params_new@other_params$new_refuge <- legacy$other_params$new_refuge
    cat("Set new_refuge flag to:", legacy$other_params$new_refuge, "\n")
}

# ========================================================================
# RUN TO STEADY STATE
# ========================================================================
# The model needs to equilibrate after reconstruction. reefSteady projects
# forward in time until biomasses stabilize.

cat("\nRunning reefSteady...\n")
params_new <- reefSteady(params_new)

# ========================================================================
# SAVE MODERNIZED VERSION
# ========================================================================
# Save both as .rds (with mizer::saveParams for full fidelity) and as .rda
# (for package data compatibility).

cat("Saving modernized params...\n")
mizer::saveParams(params_new, modern_params_rds)
save(params_new, file = modern_params_rda)
cat("Modernized params saved to:\n")
cat("  ", modern_params_rds, "\n")
cat("  ", modern_params_rda, "\n\n")

# ========================================================================
# GENERATE COMPARISON PLOTS
# ========================================================================
# Create the same plots from the modernized version so you can compare
# them side-by-side with the legacy plots to see if reconstruction worked
# or if there are systematic differences.

cat("Generating comparison plots...\n")
tryCatch({
    png(here("artifacts/modern-biomass-vs-species.png"), width = 1400, height = 900, res = 150)
    print(plotBiomassVsSpecies(params_new))
    dev.off()
}, error = function(e) cat("Could not plot modern biomass:", e$message, "\n"))

tryCatch({
    png(here("artifacts/modern-spectra-power1.png"), width = 1400, height = 900, res = 150)
    print(plotSpectra(params_new, total = TRUE, power = 1))
    dev.off()
}, error = function(e) cat("Could not plot modern spectra:", e$message, "\n"))

tryCatch({
    png(here("artifacts/modern-spectra-power2.png"), width = 1400, height = 900, res = 150)
    print(plotSpectra(params_new, total = TRUE, power = 2))
    dev.off()
}, error = function(e) cat("Could not plot modern spectra:", e$message, "\n"))

# ========================================================================
# SUMMARY
# ========================================================================
cat("\n=== RECOVERY COMPLETE ===\n")
cat("Check these files for forensic analysis:\n")
cat("  1. ", legacy_diagnostics_txt, " - Full parameter dump\n")
cat("     Look for: erepro values, reproduction_level, resource_params,\n")
cat("               refuge method, degrade flags, biomass values\n\n")
cat("  2. artifacts/legacy-plots/ - All diagnostic plots from legacy model\n")
cat("     Compare shapes/values to what you remember of 'test 2'\n\n")
cat("  3. ", legacy_dump_rds, " - Raw legacy data (species_params, interaction, etc.)\n")
cat("     Load this if you need to inspect parameters programmatically\n\n")
cat("  4. ", modern_params_rds, " - Modernized version\n")
cat("     Use this as a starting point for further calibration\n\n")
cat("Compare legacy-plots/ to modern plots to identify differences.\n")
cat("If modern plots differ significantly, the reconstruction may need adjustment.\n")
cat("Review the diagnostics file for all parameter values to identify which 'test' this was.\n")
