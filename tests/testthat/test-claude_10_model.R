test_that("claude_10_model differs from caribbean_10_model only in algae_growth and detritus external", {
    data(caribbean_10_model)
    data(claude_10_model)

    original <- caribbean_10_model
    original@other_params$algae_params$algae_growth <-
        claude_10_model@other_params$algae_params$algae_growth
    original@other_params$detritus_params$external <-
        claude_10_model@other_params$detritus_params$external
    original@time_modified <- claude_10_model@time_modified

    expect_equal(original, claude_10_model)
})

test_that("claude_10_model's algae_growth matches what tuneUR() computes fresh from caribbean_10_model", {
    data(caribbean_10_model)
    data(claude_10_model)
    expected_growth <- sum(getAlgaeConsumption(caribbean_10_model))
    expect_equal(claude_10_model@other_params$algae_params$algae_growth, expected_growth)
})

test_that("claude_10_model reaches a genuine steady state for algae and detritus, unlike caribbean_10_model", {
    data(caribbean_10_model)
    data(claude_10_model)

    for (m in list(caribbean_10_model = caribbean_10_model, claude_10_model = claude_10_model)) {
        P_A <- sum(getAlgaeProduction(m))
        c_A <- algae_consumption(m, n = m@initial_n, rates = getRates(m))
        algae_dbdt <- P_A - c_A * algae_biomass(m)

        P_D <- sum(getDetritusProduction(m))
        c_D <- detritus_consumption(m, n = m@initial_n, rates = getRates(m))
        detritus_dbdt <- P_D - c_D * detritus_biomass(m)

        if (identical(m@other_params$algae_params$algae_growth, 2000)) {
            # caribbean_10_model: known not to be tuned
            expect_gt(abs(algae_dbdt), 1)
            expect_gt(abs(detritus_dbdt), 1)
        } else {
            # claude_10_model: tuned
            expect_equal(algae_dbdt, 0)
            expect_equal(detritus_dbdt, 0, tolerance = 1e-8)
        }
    }
})

test_that("claude_10_model is a valid mizerReef object that can be projected without producing NAs", {
    data(claude_10_model)
    expect_s4_class(claude_10_model, "mizerReef")
    expect_error(mizer::validParams(claude_10_model), NA)

    sim <- project(claude_10_model, t_max = 2, progress_bar = FALSE)
    expect_false(anyNA(sim@n))
    expect_false(anyNA(sim@n_other))
})
