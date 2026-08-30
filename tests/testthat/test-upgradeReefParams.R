# upgradeReefParams() migrates params objects built with mizerReef 1.x (the
# thesis-era `main` branch: plain "MizerParams" class, direct `reef*`
# @rates_funcs overrides, flat other_params with `bad_pred` etc.) to the
# current, nested other_params / "mizerReef" S3-dispatch layout.
#
# Real mizerReef 1.x objects live on a different branch (`main`) and are not
# reachable from this package's own test suite, so these tests build a
# synthetic, structurally-accurate stand-in: start from a real object built
# with the CURRENT newReefParams() (so every core mizer slot is valid), then
# reshape its other_params/species_params/rates_funcs down into the flat
# 1.x shape documented in inst/to-do-list.txt (verified there against the
# real bundled `main`-branch karpata_model.rda/bonaire_model.rda objects).
make_fake_v1 <- function(method = "binned",
                         method_params = data.frame(
                             start_L = c(0, 20), end_L = c(20, 60),
                             prop_protect = c(0.9, 0.5)
                         ),
                         degrade = FALSE) {
    sp <- data.frame(
        species = c("a", "b", "c"),
        w_max = c(50, 500, 5000),
        w_mat = c(5, 50, 500),
        k_vb = c(1, 1, 1),
        a = c(0.01, 0.01, 0.01), b = c(3, 3, 3),
        interaction_algae = c(0.5, 0, 0),
        interaction_detritus = c(0, 0.5, 0),
        refuge_user = c(TRUE, TRUE, FALSE),
        blocked_pred = c(FALSE, FALSE, TRUE),
        satiation = c(TRUE, FALSE, FALSE)
    )
    p <- suppressWarnings(newReefParams(sp,
        method = method, method_params = method_params,
        degrade = degrade,
        trajectory = if (degrade) "bleaching" else NULL,
        bleach_time = 3,
        deg_scale = if (degrade) matrix(c(1, 0.9, 0.8, 0.7), nrow = 1) else NULL
    ))

    op <- p@other_params
    new_op <- list(
        include_sen_mort = op$include_sen_mort,
        refuge_params = data.frame(
            method = op$refuge_params$method,
            a_bar = op$refuge_params$a_bar, b_bar = op$refuge_params$b_bar,
            w_settle = op$refuge_params$w_settle,
            max_protect = op$refuge_params$max_protect,
            tau = op$refuge_params$tau
        ),
        method_params = op$refuge_params$method_params,
        refuge = op$refuge_params$refuge,
        bin.id = op$refuge_params$bin.id,
        refuge_lengths = op$refuge_params$refuge_lengths,
        carry_capacity = FALSE,
        initial_algae_growth = 2000,
        algae_capacity = op$algae$capacity,
        sen_decomp = op$detritus$sen_decomp,
        ext_decomp = op$detritus$ext_decomp,
        initial_d_external = 1,
        detritus_capacity = op$detritus$capacity,
        ext_mort_params = op$ext_mort_params,
        algae = list(rho = op$algae$rho, growth = op$algae$growth),
        detritus = list(
            rho = op$detritus$rho, sen_decomp = op$detritus$sen_decomp,
            ext_decomp = op$detritus$ext_decomp, external = op$detritus$external
        ),
        new_refuge = FALSE,
        degrade = degrade
    )
    if (degrade) {
        new_op$t_bleach <- op$refuge_params$t_bleach
        new_op$trajectory <- op$refuge_params$trajectory
        new_op$deg_scale <- op$refuge_params$deg_scale
    }
    p@other_params <- new_op
    names(p@species_params)[names(p@species_params) == "blocked_pred"] <- "bad_pred"
    p@rates_funcs$Rates <- "reefRates"
    p@rates_funcs$Encounter <- "reefEncounter"
    p@rates_funcs$FeedingLevel <- "reefFeedingLevel"
    p@rates_funcs$PredMort <- "reefPredMort"
    p@rates_funcs$Mort <- "reefMort"
    p@extensions <- character(0)
    class(p) <- "MizerParams"
    p
}

# upgradeReefParams() emits two informative warning()s per migration; nesting
# expect_warning() calls to check both is fragile (the inner call consumes/
# muffles its match before the outer one runs), so collect every warning
# message up front instead and check the return value and the messages
# separately.
capture_result_and_warnings <- function(expr) {
    messages <- character(0)
    result <- withCallingHandlers(
        expr,
        warning = function(w) {
            messages <<- c(messages, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
    )
    list(result = result, warnings = messages)
}

test_that("upgradeReefParams migrates a fake 1.x object to a valid mizerReef object", {
    fake <- make_fake_v1()
    out <- capture_result_and_warnings(upgradeReefParams(fake))

    expect_s3_class(out$result, "mizerReef")
    expect_true(inherits(mizer::validParams(out$result), "mizerReef"))
    expect_match(out$warnings, "rate-function overrides", all = FALSE)
    expect_match(out$warnings, "not guaranteed to be at a genuine joint steady state",
                 all = FALSE)
})

test_that("upgradeReefParams renames bad_pred to blocked_pred", {
    fake <- make_fake_v1()
    result <- suppressWarnings(upgradeReefParams(fake))
    expect_true("blocked_pred" %in% names(result@species_params))
    expect_false("bad_pred" %in% names(result@species_params))
    expect_equal(result@species_params$blocked_pred, fake@species_params$bad_pred)
})

test_that("upgradeReefParams drops bad_pred without clobbering an existing blocked_pred", {
    fake <- make_fake_v1()
    fake@species_params$blocked_pred <- rep(TRUE, nrow(fake@species_params))
    result <- suppressWarnings(upgradeReefParams(fake))
    expect_false("bad_pred" %in% names(result@species_params))
    expect_true(all(result@species_params$blocked_pred))
})

test_that("upgradeReefParams resets the old direct rate-function overrides", {
    fake <- make_fake_v1()
    result <- suppressWarnings(upgradeReefParams(fake))
    expect_equal(result@rates_funcs$Rates, "mizerRates")
    expect_equal(result@rates_funcs$Encounter, "mizerEncounter")
    expect_equal(result@rates_funcs$FeedingLevel, "mizerFeedingLevel")
    expect_equal(result@rates_funcs$PredMort, "mizerPredMort")
    expect_equal(result@rates_funcs$Mort, "mizerMort")
})

test_that("upgradeReefParams preserves tuned algae/detritus values instead of resetting them", {
    fake <- make_fake_v1()
    tuned_growth <- fake@other_params$algae$growth
    tuned_external <- fake@other_params$detritus$external
    result <- suppressWarnings(upgradeReefParams(fake))

    expect_equal(result@other_params$algae$growth, tuned_growth)
    expect_equal(result@other_params$detritus$external, tuned_external)
})

test_that("upgradeReefParams builds the expected algae/detritus/refuge schema", {
    fake <- make_fake_v1()
    result <- suppressWarnings(upgradeReefParams(fake))

    expect_setequal(names(result@other_params$algae), c("rho", "capacity", "growth"))
    expect_null(result@other_params$algae_params)
    expect_setequal(names(result@other_params$detritus),
                    c("rho", "capacity", "sen_decomp", "ext_decomp", "external"))
    expect_null(result@other_params$detritus_params)
    expect_true(is.list(result@other_params$refuge_params))
    expect_false(is.data.frame(result@other_params$refuge_params))
    expect_true(result@other_params$refuge_params$use_dummy_fish_bins)
})

test_that("upgradeReefParams removes the old flat top-level fields", {
    fake <- make_fake_v1()
    result <- suppressWarnings(upgradeReefParams(fake))

    stale <- c(
        "method_params", "refuge", "bin.id", "refuge_lengths",
        "carry_capacity", "algae_capacity", "detritus_capacity",
        "sen_decomp", "ext_decomp", "initial_algae_growth",
        "initial_d_external"
    )
    expect_false(any(stale %in% names(result@other_params)))
})

test_that("upgradeReefParams leaves core mizer slots untouched", {
    fake <- make_fake_v1()
    result <- suppressWarnings(upgradeReefParams(fake))

    core_slots <- c("w", "dw", "w_full", "initial_n", "search_vol", "psi",
                    "intake_max", "metab", "interaction")
    for (s in core_slots) {
        expect_equal(result[[s]], fake[[s]], info = s)
    }
})

test_that("upgradeReefParams produces a params object that projects without error", {
    fake <- make_fake_v1()
    result <- suppressWarnings(upgradeReefParams(fake))
    expect_no_error(project(result, t_max = 1, progress_bar = FALSE))
})

test_that("upgradeReefParams works for each refuge method", {
    method_params_by_method <- list(
        sigmoidal = data.frame(L_refuge = 20, prop_protect = 0.9),
        binned = data.frame(start_L = c(0, 20), end_L = c(20, 60), prop_protect = c(0.9, 0.5)),
        competitive = data.frame(start_L = c(0, 20), end_L = c(20, 60), refuge_density = c(5, 2)),
        noncomplex = NULL
    )
    for (method in names(method_params_by_method)) {
        fake <- make_fake_v1(method = method, method_params = method_params_by_method[[method]])
        result <- suppressWarnings(upgradeReefParams(fake))
        expect_equal(result@other_params$refuge_params$method, method, info = method)
        expect_s3_class(result, "mizerReef")
    }
})

test_that("upgradeReefParams relocates degradation trajectory settings and warns", {
    fake <- make_fake_v1(degrade = TRUE)
    out <- capture_result_and_warnings(upgradeReefParams(fake))

    expect_match(out$warnings, "rate-function overrides", all = FALSE)
    expect_match(out$warnings, "habitat degradation enabled", all = FALSE)
    expect_match(out$warnings, "not guaranteed to be at a genuine joint steady state",
                 all = FALSE)

    result <- out$result
    expect_true(result@other_params$refuge_params$degrade)
    expect_equal(result@other_params$refuge_params$t_bleach, fake@other_params$t_bleach)
    expect_equal(result@other_params$refuge_params$trajectory, fake@other_params$trajectory)
    expect_equal(result@other_params$refuge_params$deg_scale, fake@other_params$deg_scale)
})

test_that("upgradeReefParams leaves an already-current mizerReef object unchanged", {
    data(caribbean_3_model)
    messages <- character(0)
    result <- withCallingHandlers(
        upgradeReefParams(caribbean_3_model),
        message = function(m) {
            messages <<- c(messages, conditionMessage(m))
            invokeRestart("muffleMessage")
        }
    )
    expect_match(messages, "already a mizerReef object", all = FALSE)
    expect_s3_class(result, "mizerReef")
    expect_equal(result@other_params, caribbean_3_model@other_params)
})

test_that("upgradeReefParams errors on an object that isn't mizerReef-shaped at all", {
    data(NS_species_params)
    plain <- suppressWarnings(suppressMessages(
        newMultispeciesParams(NS_species_params, no_w = 30)
    ))
    expect_error(upgradeReefParams(plain), "does not look like a mizerReef 1.x")
})
