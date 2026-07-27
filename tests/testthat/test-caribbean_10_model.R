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
