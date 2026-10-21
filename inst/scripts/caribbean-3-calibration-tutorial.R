# Calibrating the 3-group Caribbean reef model from literature values
# A worked example to follow when calibrating your own mizerReef model
#
# This script starts from literature parameter values and repeats the
# calibration of caribbean_3_model with the standard package steps from
# vignette("steady-state-recipe"), plus a small solver (Step 8) for the
# targets that those steps cannot reach one at a time. It ends very close to
# the shipped model. The shipped model itself is built without any
# calibration by caribbean-3-calibration-clean.r, from parameter tables that
# already hold the calibrated values.
#
# The model follows the structure of Rogers et al. (2018): herbivores eat only
# algae, invertebrates eat only detritus, and only predators eat plankton, up
# to 10^-1.5 g (the smallest fish size). Fish compete for refuge holes
# (Karpata reef, Bonaire) once they are 2.5 cm long.
#
# Targets at steady state:
#   - biomass above 3.125 g of 107 (predators), 34 (herbivores) and
#     40 (invertebrates) g/m^2, from the FORCE surveys in Bonaire
#   - about 10 herbivore juveniles (0.03-1.6 g, or 1.1-4 cm) per m^2: Rogers
#     et al. (2018) set recruitment to give about 10 small fish per m^2 in
#     each fish group
#   - adult predators eat 33% predators, 17% herbivores and 50%
#     invertebrates: Rogers' fractions of feeding time, used here as diet
#     shares
#   - predators mature at age 2 (age_mat in the species table)
#   - algae standing biomass of 350 g/m^2 and a detritus lifetime of
#     3.85 hours
#
# Run it from the package source folder. It takes a few minutes; the solver
# in Steps 8 and 9 is the slow part.

# Load required libraries
library(mizer)
library(mizerReef)
library(here)

# Step 0: targets and a function to check them -------------------------------
targets <- c(B_herb = 34, B_pred = 107, B_inv = 40, J_herb = 10,
             diet_pred = 0.33, diet_herb = 0.17, age_mat_pred = 2)

# Adult predators' diet: the share of predators, herbivores and invertebrates
# in what predators larger than their maturity size eat (plankton left out)
adult_predator_diet <- function(params) {
  diet <- getDiet(params, proportion = FALSE)
  adult <- params@w >= params@species_params["predators", "w_mat"]
  eaten <- colSums(diet["predators", adult, ] *
                     params@initial_n["predators", adult] * params@dw[adult])
  prey <- eaten[c("predators", "herbivores", "inverts")]
  prey / sum(prey)
}

# The model's value for each target
calibration_outputs <- function(params) {
  biomass <- getBiomass(params, use_cutoff = TRUE)
  diet <- adult_predator_diet(params)
  c(B_herb = biomass[["herbivores"]],
    B_pred = biomass[["predators"]],
    B_inv = biomass[["inverts"]],
    J_herb = getN(params, min_w = 10^-1.5, max_w = 1.6)[["herbivores"]],
    diet_pred = diet[["predators"]],
    diet_herb = diet[["herbivores"]],
    age_mat_pred = age_mat(params)[["predators"]])
}
check_targets <- function(params) {
  print(signif(rbind(model = calibration_outputs(params)[names(targets)],
                     target = targets), 3))
}

# Step 1: load the parameter tables and go back to literature values ---------
caribbean_3_species <- read.csv(here("inst/data-csv/caribbean_3_species.csv"))
caribbean_3_interaction <- read.csv(
  here("inst/data-csv/caribbean_3_interaction.csv"), row.names = 1)
karpata_refuge <- read.csv(here("inst/data-csv/karpata_refuge.csv"))
tuning_profile <- read.csv(here("inst/data-csv/tuning_profile.csv"))

# These tables hold the calibrated values that the shipped model is built
# from. To calibrate from scratch, take the calibrated columns out again so
# that newReefParams() works out its own defaults (gamma, h, erepro, R_max),
# and put back the starting values used before calibration:
species <- caribbean_3_species
species[c("gamma", "h", "erepro", "R_max")] <- NULL
# Standard metabolism: the values in the species table before calibration
# (no source recorded; they equal Rogers' growth efficiencies, 0.15 and 0.1,
# but ks is a metabolic rate, not an efficiency)
species$ks <- c(0.15, 0.1, 0.1)
# Herbivores eat algae and invertebrates eat detritus at full strength; these
# only set the scale of algae and detritus, which Step 10 sets properly
species$interaction_algae <- c(0, 1, 0)
species$interaction_detritus <- c(0, 0, 1)
# Predators' diet weights: Rogers et al. (2018) Table S2
interaction <- caribbean_3_interaction
interaction["predators", ] <- c(0.33, 0.17, 0.5)
# Everything else stays as in the table. From Rogers et al. (2018) Table S2:
# fish recruit at w_min = 10^-1.5 g, search volume scales with w^0.75
# (q = 0.75), only predators eat plankton (interaction_resource 1/0/0).
# Invertebrates keep mizer's default w_min of 0.001 g (Rogers: 10^-4 g).
species[, c("species", "w_min", "q", "ks", "interaction_resource",
            "interaction_algae", "interaction_detritus", "satiation")]

# Step 2: set up the model with the simple (binned) refuge -------------------
# The competitive refuge depends on fish numbers, which makes it hard to find
# a first steady state, so calibration starts with a constant refuge profile
# (60% of fish protected at every size). See vignette("steady-state-recipe").
#
# Plankton: Rogers' plankton only reaches 10^-1.5 g, too small to feed
# predators larger than a few grams. With mizer's starting abundances fish
# are far too scarce to feed them instead, so predators would stop growing
# at 30-100 g and never reach maturity. We start with
# plankton up to 1 g and bring it down to 10^-1.5 g in Step 6, once
# predators have learnt to eat fish.
params <- newReefParams(
  species_params = species,
  interaction = interaction,
  method = "binned",
  method_params = tuning_profile,
  # Fish use refuge from 2.5 cm (0.39 g); this setting carries over to the
  # competitive refuge in Step 7
  w_settle = 0.025 * 2.5^3,
  w_pp_cutoff = 1)

# newReefParams() warns that kappa is still at mizer's placeholder value
# (1e11). Herbivores here eat only algae, and algae production is a fixed
# 2000 g/m^2/year, so the fish have to be near their real abundance before
# the first steady state: at the placeholder scale there are about 1e8 g/m^2
# of herbivores sharing that algae, and they starve. So calibrate the overall
# scale straight away.
params <- calibrateReefBiomass(params)

# Moderate density dependence in reproduction helps the model settle
params <- setBevertonHolt(params, reproduction_level = 0.5)

# Step 3: find a first steady state ------------------------------------------
params <- params |>
  reefSteady() |> reefSteady() |> reefSteady() |>
  reefSteady() |> reefSteady() |> reefSteady()

# Step 4: match biomasses and age at maturity --------------------------------
# One round of matching. matchBiomasses() and matchReefGrowth() keep erepro
# fixed, so the Beverton-Holt reproduction level drifts, sometimes all the
# way to 1, where reefSteady() stops with an error. Setting it back to 0.5
# each round prevents that; it changes erepro and R_max but not the
# abundances, because reefSteady() holds recruitment fixed anyway.
match_round <- function(params) {
  params |>
    matchBiomasses() |> matchReefGrowth() |>
    setBevertonHolt(reproduction_level = 0.5) |>
    reefSteady()
}
params <- calibrateReefBiomass(params) |> reefSteady()
n_iter <- 10
for (i in 1:n_iter) {
  params <- match_round(params)
}
plotBiomassObservedVsModel(params)
age_mat(params) # predators and herbivores should be at 2
check_targets(params)

# Step 5: shift predators from plankton to fish ------------------------------
# The biomasses match, but adult predators eat almost nothing but plankton:
plotDiet(params, species = "predators")
# scaleDownPlankton() divides plankton by a factor and multiplies every
# search volume by it, so fish prey become that many times more important.
# One big step throws the model too far from steady state, so take five
# steps of 4 (about 1000 in total) and re-match after each.
for (i in 1:5) {
  params <- scaleDownPlankton(params, factor = 4) |> reefSteady()
  for (j in 1:4) {
    params <- match_round(params)
  }
}
plotDiet(params, species = "predators") # now mostly fish

# Step 6: plankton only up to 10^-1.5 g, as in Rogers et al. (2018) ----------
# Now that predators eat fish, the plankton ceiling can come down in one go;
# it takes about 15 rounds of re-matching to settle again
resource_params(params)$w_pp_cutoff <- 10^-1.5
# remove the plankton above the new ceiling
initialNResource(params) <- initialNResource(params) *
  (resource_capacity(params) > 0)
params <- reefSteady(params)
for (i in 1:15) {
  params <- match_round(params)
}
check_targets(params)

# Step 7: switch to the competitive refuge -----------------------------------
# Now fish compete for the holes counted on Karpata reef. setRefuge() is used
# rather than newRefuge(), because newRefuge() would put w_settle back to its
# default of 0.1 g.
params <- setRefuge(params, method = "competitive",
                    method_params = karpata_refuge,
                    w_settle = 0.025 * 2.5^3) |>
  getRefuge()
# The refuge changes how much energy goes into eggs, and with it the
# reproduction level; set it back to 0.5
params <- setBevertonHolt(params, reproduction_level = 0.5)
# Let the model settle with recruitment held fixed. Don't use
# matchBiomasses() from here on: it scales every size of a group by the same
# factor, juveniles included. Here that leaves almost no juveniles (about
# 0.004 per m^2 in earlier tests), the holes stop being contested, and the
# solver below could not recover from that state.
params <- params |> reefSteady() |> reefSteady() |> reefSteady()
plotRefugeProfile(params)
check_targets(params) # herbivores are now far too many

# Step 8: solve for the coupled targets --------------------------------------
# The targets now depend on each other: predators eat herbivores, herbivores
# and young predators share the holes, and predators eat invertebrates. Fixing
# one target at a time undoes the others, so a solver changes seven knobs
# together to hit the seven targets:
#   - predators' interaction with herbivores, predators and invertebrates
#   - the recruitment of each group (recruits per m^2 per year)
#   - the plankton abundance (kappa)
# After every change reefSteady() settles the model with recruitment held.
# Predators' search volume (gamma) is not a knob: gamma times the
# interaction is what matters, so gamma and the three interactions together
# would be one knob too many.

# Recruitment at an unchanged Beverton-Holt reproduction level: erepro and
# R_max are set so that recruitment is exactly `recruitment`
set_recruitment <- function(params, recruitment, level) {
  rdi <- getRDI(params)
  sp <- species_params(params)
  for (s in names(recruitment)) {
    sp[s, "erepro"] <- sp[s, "erepro"] * recruitment[[s]] /
      (1 - level[[s]]) / rdi[[s]]
    sp[s, "R_max"] <- recruitment[[s]] / level[[s]]
  }
  species_params(params, recalculate = FALSE) <- sp
  params
}
get_knobs <- function(params) {
  recruitment <- getRDD(params)
  c(pred_on_herb = params@interaction["predators", "herbivores"],
    pred_on_pred = params@interaction["predators", "predators"],
    pred_on_inv = params@interaction["predators", "inverts"],
    R_herb = recruitment[["herbivores"]],
    R_inv = recruitment[["inverts"]],
    R_pred = recruitment[["predators"]],
    kappa = params@resource_params$kappa)
}
set_knobs <- function(params, knobs) {
  level <- getReproductionLevel(params)
  inter <- params@interaction
  inter["predators", "herbivores"] <- knobs[["pred_on_herb"]]
  inter["predators", "predators"] <- knobs[["pred_on_pred"]]
  inter["predators", "inverts"] <- knobs[["pred_on_inv"]]
  params <- setInteraction(params, interaction = inter)
  plankton_factor <- knobs[["kappa"]] / params@resource_params$kappa
  resource_params(params)$kappa <- knobs[["kappa"]]
  initialNResource(params) <- initialNResource(params) * plankton_factor
  set_recruitment(params, level = level,
                  recruitment = c(herbivores = knobs[["R_herb"]],
                                  inverts = knobs[["R_inv"]],
                                  predators = knobs[["R_pred"]]))
}
settle <- function(params) {
  params |>
    reefSteady(progress_bar = FALSE, info_level = 0) |>
    reefSteady(progress_bar = FALSE, info_level = 0) |>
    reefSteady(progress_bar = FALSE, info_level = 0)
}

# A damped Newton solver. Each iteration runs the model once per knob, with
# that knob changed by 10%, to see how the targets respond (the Jacobian).
# It then moves all knobs part of the way (damp) towards where every target
# would be met, by at most a factor `cap` in one step. `limits` stops chosen
# knobs from moving more than that factor away from their starting values.
# Errors are on a log scale: 0.01 means 1%.
solve_targets <- function(params, targets, limits = NULL, n_iter = 60,
                          tol = 0.01, damp = 0.5, cap = 2, h = 0.1,
                          cores = if (.Platform$OS.type == "windows") 1 else 4) {
  log_error <- function(p) log(calibration_outputs(p)[names(targets)] / targets)
  x <- log(get_knobs(params))
  lower <- x[names(limits)] - log(limits)
  upper <- x[names(limits)] + log(limits)
  err <- log_error(params)
  if (max(abs(err)) > log(10)) {
    message("Starting far from the targets (", names(which.max(abs(err))),
            " is off by a factor ", signif(exp(max(abs(err))), 2), "). ",
            "The solver moves at most a factor ", cap, " per step, so it may ",
            "take many steps. If it fails, get closer with the standard ",
            "steps first.")
  }
  history <- NULL
  best <- params
  best_err <- err
  for (it in 0:n_iter) {
    history <- rbind(history, c(iteration = it, exp(x),
                                calibration_outputs(params)[names(targets)],
                                worst_error = max(abs(err))))
    message("Solver iteration ", it, ": worst target off by ",
            signif(100 * max(abs(err)), 2), "%")
    if (max(abs(err)) < max(abs(best_err))) {
      best <- params
      best_err <- err
    }
    if (max(abs(err)) < tol || it == n_iter) break
    # How each target responds to each knob
    columns <- parallel::mclapply(names(x), function(k) {
      x_k <- x
      x_k[[k]] <- x_k[[k]] + h
      p_k <- tryCatch(settle(set_knobs(params, exp(x_k))),
                      error = function(e) NULL)
      if (is.null(p_k)) return(rep(NA, length(err)))
      (log_error(p_k) - err) / h
    }, mc.cores = cores)
    J <- do.call(cbind, columns)
    dimnames(J) <- list(names(targets), names(x))
    if (any(!is.finite(J))) {
      message("Stopped: a model run failed while testing the knobs.")
      break
    }
    # If two knobs change the targets in the same way, the solver can't
    # tell them apart (the Jacobian is singular)
    s <- svd(J)
    if (min(s$d) < 1e-6 * max(s$d)) {
      same <- names(x)[abs(s$v[, which.min(s$d)]) > 0.3]
      message("Stopped: the knobs ", paste(same, collapse = ", "),
              " change the targets in the same way. Either they are ",
              "redundant (drop one, or add a target that separates them), ",
              "or the model has moved to a state where some targets no ",
              "longer respond (see the solver history).")
      break
    }
    step <- -damp * solve(J, err)
    step <- step * min(1, log(cap) / max(abs(step)))
    # Keep limited knobs within their limits
    x_new <- x + step
    x_new[names(limits)] <- pmin(pmax(x_new[names(limits)], lower), upper)
    step <- x_new - x
    # Take the step; if the model fails, try half the step (up to 3 times)
    p_new <- NULL
    for (halving in 0:3) {
      p_new <- tryCatch(settle(set_knobs(params, exp(x + step))),
                        error = function(e) NULL)
      if (!is.null(p_new) && all(is.finite(log_error(p_new)))) break
      p_new <- NULL
      step <- step / 2
    }
    if (is.null(p_new)) {
      message("Stopped: the model failed even with a small step.")
      break
    }
    params <- p_new
    x <- x + step
    err <- log_error(params)
  }
  missed <- names(targets)[abs(best_err) >= tol]
  if (length(missed) > 0) {
    out <- calibration_outputs(best)
    message("Not every target was met. The closest state found misses: ",
            paste0(missed, " (", signif(out[missed], 3), " vs ",
                   targets[missed], ")", collapse = ", "))
  }
  attr(best, "solver_history") <- history
  best
}

# Limits: predators' interactions may move at most a factor 10 from Rogers'
# values, and plankton at most a factor 10 from where Steps 5 and 6 left it.
# Without them the solver can end up somewhere else:
#   - a dead end where predators eat almost only their own young
#     (interaction with predators above 10, thousands of predator recruits
#     per m^2) and the targets stop responding;
#   - a second solution that meets every target, but in which plankton is
#     about 100 times more abundant, young predators live on it until they
#     are about 100 g, and there are almost no small predators (0.2 recruits
#     per m^2 per year, against about 10 small fish per m^2 per group in
#     Rogers et al. 2018).
limits <- c(pred_on_herb = 10, pred_on_pred = 10, pred_on_inv = 10,
            kappa = 10)
params <- solve_targets(params, targets, limits = limits)
signif(attr(params, "solver_history"), 3)
check_targets(params)

# Step 9: give algae a carrying capacity -------------------------------------
# setURcapacity() sets the capacity to 1.5 times the current biomass, for
# algae and detritus. Detritus is not a growing population, so it gets no
# capacity (Inf). Once herbivores graze it, algae settles at about 60% of
# its earlier biomass (so the capacity ends up about 2.5 times the algae
# biomass) and herbivores get less food. Solve again; this time the solver
# starts close to the targets.
params <- setURcapacity(params, cap = 1.5)
params@other_params$detritus$capacity <- Inf
params <- params |> reefSteady() |> reefSteady() |> reefSteady()
check_targets(params)
params <- solve_targets(params, targets, limits = limits)
check_targets(params)

# Step 10: put algae and detritus on their literature scale ------------------
# Up to now algae and detritus biomasses have been on whatever scale the
# rescaling steps left them (tiny). rescale_algae() and detritus_lifetime<-
# multiply a biomass by a factor and divide the consumers' encounter rate by
# it, so the fish see exactly the same food. The algae carrying capacity has
# to move with the algae biomass (rescale_algae() does not do this itself).
# Targets: algae 350 g/m^2, the midpoint of the 150-550 g wet weight/m^2
# range used for the 10-group model (converted from a Caribbean estimate of
# about 10-100 g dry weight/m^2); detritus lifetime 3.85 hours, the midpoint
# of the 2.8-4.9 hour residence time measured at Palmyra Atoll by Max et al.
# (2013).
algae_factor <- 350 / algae_biomass(params)
params <- rescale_algae(params, algae_factor)
params@other_params$algae$capacity <-
  params@other_params$algae$capacity * algae_factor
detritus_lifetime(params) <- 3.85 / (24 * 365)

# Step 11: final checks ------------------------------------------------------
# Is it a stable steady state? Project 100 years with Beverton-Holt
# reproduction live: the biomasses should barely move (reefSteady() stops
# within about 1% of the steady state). Keep the end state, as
# caribbean-3-calibration-clean.r does.
sim <- project(params, t_max = 100, t_save = 5, progress_bar = FALSE)
plotBiomass(sim)
biomass <- getBiomass(sim)[, params@species_params$species]
max(abs(biomass[nrow(biomass), ] / biomass[1, ] - 1)) # close to 0
params <- finalParams(sim)

check_targets(params)
plotBiomassObservedVsModel(params)
plotDiet(params, species = "predators")
age_mat(params)
algae_biomass(params)                  # target 350 g/m^2
params@other_params$algae$capacity     # algae carrying capacity (g/m^2)
detritus_lifetime(params) * 24 * 365   # target 3.85 hours
# erepro above 1 is unrealistic: the species would need to turn more than
# all of its reproductive energy into eggs. A lower reproduction level (more
# density dependence) gives a lower erepro and the same steady state:
# setBevertonHolt(params, reproduction_level = ...)
params@species_params$erepro
getReproductionLevel(params)

# Step 12: compare with the shipped model ------------------------------------
shipped <- mizerReef::caribbean_3_model
signif(rbind(this_script = calibration_outputs(params),
             shipped = calibration_outputs(shipped)), 3)
# Predators' encounter with each prey depends on gamma times the interaction,
# so compare that product rather than either on its own
key_values <- function(p) {
  sp <- p@species_params
  c(sp["predators", "gamma"] * p@interaction["predators", ],
    kappa = p@resource_params$kappa,
    setNames(getRDD(p), paste0("R_", sp$species)),
    setNames(sp$h, paste0("h_", sp$species)),
    setNames(sp$ks, paste0("ks_", sp$species)),
    setNames(sp$erepro, paste0("erepro_", sp$species)),
    algae_capacity = p@other_params$algae$capacity,
    detritus_external = p@other_params$detritus$external,
    age_mat_herb = age_mat(p)[["herbivores"]])
}
signif(rbind(this_script = key_values(params),
             shipped = key_values(shipped)), 3)

# Step 13: save --------------------------------------------------------------
# Saved to a temporary file; change the path to keep it. The shipped
# caribbean_3_model is built by caribbean-3-calibration-clean.r instead.
saveParams(params, file.path(tempdir(), "caribbean_3_calibrated.rds"))
