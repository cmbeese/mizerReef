test_that("caribbean_10_model is a valid mizerReef object", {
    data(caribbean_10_model)
    params <- caribbean_10_model
    expect_s4_class(params, "mizerReef")
    expect_error(mizer::validParams(params), NA)
    expect_equal(nrow(params@species_params), 10)
})

test_that("caribbean_10_model's rates_funcs hold mizer's default names, not stale reef-specific strings", {
    # Regression check for the "stale bundled .rda metadata from an old
    # porting event" bug class fixed in a prior session (rates_funcs used
    # to hold literal "reefEncounter"/"reefFeedingLevel"/"reefPredMort"/
    # "reefMort" strings, which broke composability with other stacked
    # mizer extensions). The reef-specific behaviour is applied via S3
    # dispatch on the object's class instead (see the next test).
    data(caribbean_10_model)
    rf <- caribbean_10_model@rates_funcs
    expect_equal(rf$Encounter, "mizerEncounter")
    expect_equal(rf$FeedingLevel, "mizerFeedingLevel")
    expect_equal(rf$PredMort, "mizerPredMort")
    expect_equal(rf$Mort, "mizerMort")
})

test_that("caribbean_10_model's refuge-blocking is actually applied via S3 dispatch", {
    data(caribbean_10_model)
    params <- caribbean_10_model
    n <- params@initial_n
    n_pp <- params@initial_n_pp
    n_other <- params@initial_n_other

    via_getEncounter <- unclass(mizer::getEncounter(params))
    via_direct_reef <- reefEncounter(params, n, n_pp, n_other, t = 0)
    expect_equal(via_getEncounter, via_direct_reef, ignore_attr = TRUE)

    via_bare_mizer <- mizer::mizerEncounter(params, n, n_pp, n_other, t = 0)
    expect_false(isTRUE(all.equal(unname(via_getEncounter), unname(via_bare_mizer))))
})

test_that("caribbean_10_model has use_dummy_fish_bins set explicitly", {
    # Regression check for the use_dummy_fish_bins default-mismatch bug
    # that made caribbean_3_model's predation refuge completely
    # non-functional -- confirm caribbean_10_model doesn't rely on the
    # (now-fixed) default either.
    data(caribbean_10_model)
    expect_true(isTRUE(caribbean_10_model@other_params$refuge_params$use_dummy_fish_bins))
})

test_that("caribbean_10_model's method_params matches the bundled karpata_refuge data object", {
    data(caribbean_10_model)
    data(karpata_refuge)
    expect_equal(
        caribbean_10_model@other_params$refuge_params$method_params,
        karpata_refuge,
        ignore_attr = TRUE
    )
})

test_that("caribbean_10_model has no missing or negative abundances/parameters", {
    data(caribbean_10_model)
    params <- caribbean_10_model
    expect_false(anyNA(params@initial_n))
    expect_false(anyNA(params@initial_n_pp))
    expect_true(all(params@initial_n >= 0))
    expect_false(anyNA(params@species_params$a))
    expect_false(anyNA(params@species_params$b))
    expect_false(anyNA(params@species_params$refuge_user))
})

test_that("caribbean_10_model can be projected forward without producing NAs", {
    data(caribbean_10_model)
    sim <- project(caribbean_10_model, t_max = 2, progress_bar = FALSE)
    expect_false(anyNA(sim@n))
    expect_false(anyNA(sim@n_pp))
})

test_that("caribbean_10_model survives a multi-year live fishing simulation without negative resource biomass", {
    # Regression check for Finding 1: caribbean_10_model's fixed, large
    # negative detritus external flux used to drive total detritus
    # production negative within a fraction of a year once effort = 1
    # fishing started shifting fish abundance away from the tuned steady
    # state, eventually producing negative and then NaN detritus biomass
    # that crashed project() outright (see detritus_dynamics()'s Details
    # for the mechanism, now fixed by flooring production at zero). This is
    # the exact reproduction scenario from that diagnosis, checked at fine
    # time resolution since annual snapshots alone missed the excursion the
    # first time.
    data(caribbean_10_model)
    sim <- project(caribbean_10_model,
        effort = 1, t_max = 20, t_save = 0.1, progress_bar = FALSE
    )
    n_other <- mizer::NOther(sim)
    expect_false(anyNA(n_other))
    expect_true(all(unlist(n_other) >= 0))
})

test_that("caribbean_10_model's algae and detritus are genuinely at steady state", {
    # Regression check: an earlier version of this bundled object (now
    # archived at inst/archive/caribbean_10_model_untuned.rda) had
    # algae_growth left at newReefParams()'s raw untuned default (2000),
    # so algae/detritus were confirmed NOT at steady state (dB/dt far from
    # zero). The current object is tuneUR()'s output, which should give
    # dB/dt = 0 for both resources.
    data(caribbean_10_model)
    m <- caribbean_10_model

    P_A <- sum(getAlgaeProduction(m))
    c_A <- algae_consumption(m, n = m@initial_n, rates = getRates(m))
    expect_equal(P_A - c_A * algae_biomass(m), 0)

    P_D <- sum(getDetritusProduction(m))
    c_D <- detritus_consumption(m, n = m@initial_n, rates = getRates(m))
    expect_equal(P_D - c_D * detritus_biomass(m), 0, tolerance = 1e-8)
})

test_that("caribbean_10_model's detritus external flux is stored under the correct field name", {
    # Guard against the same stale-field-name class of bug as Finding 2
    # (caribbean_3_model's d_external): caribbean_10_model briefly carried a
    # dead d_external duplicate alongside the correctly-named external field,
    # left over from before it was tuneUR()'d. Only external should remain.
    data(caribbean_10_model)
    dp_names <- names(caribbean_10_model@other_params$detritus)
    expect_true("external" %in% dp_names)
    expect_false("d_external" %in% dp_names)
    expect_false(is.null(caribbean_10_model@other_params$detritus$external))
})

test_that("caribbean_10_model's algae/detritus consumption matrix is stored under the correct field name", {
    # Same migration-gap bug class as Finding 2, but for the algae/detritus
    # consumption matrix: newReefParams() used to write it to other_params$
    # algae_params$rho_algae / detritus_params$rho_detritus, while every
    # consumer reads/writes the bare `rho` field, only working by accident
    # of R's `$` partial name-matching. Guard the exact field name.
    data(caribbean_10_model)
    ap_names <- names(caribbean_10_model@other_params$algae)
    dp_names <- names(caribbean_10_model@other_params$detritus)
    expect_true("rho" %in% ap_names)
    expect_false("rho_algae" %in% ap_names)
    expect_true("rho" %in% dp_names)
    expect_false("rho_detritus" %in% dp_names)
})

test_that("caribbean_10_model no longer carries a duplicate algae_params/detritus_params structure", {
    # Regression check for the fix consolidating other_params$algae_params/
    # detritus_params onto the mizer-canonical other_params$algae/detritus:
    # only one structure should exist per resource.
    data(caribbean_10_model)
    expect_null(caribbean_10_model@other_params$algae_params)
    expect_null(caribbean_10_model@other_params$detritus_params)
})
