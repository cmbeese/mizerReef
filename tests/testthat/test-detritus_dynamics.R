test_that("detritus_dynamics runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(detritus_dynamics(params,
        n = params@initial_n,
        n_other = params@initial_n_other,
        rates = getRates(params),
        dt = 1), NA)
})

test_that("detritus_dynamics floors negative production at zero", {
    # Regression check for Finding 1: a fixed, large negative `external` flux
    # (as tuneUR() can produce -- see test-tuneUR.R's "negative external
    # detritus flux" test) can drive total production (feces + decomp +
    # external) negative once fish abundance shifts during a live simulation,
    # which previously broke the analytic update's non-negativity guarantee
    # and eventually produced NaN detritus biomass (see detritus_dynamics()'s
    # Details for the full mechanism). Production should be floored at zero
    # before entering the update, so the result should match the
    # zero-production formula exactly, not go negative.
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$detritus_params$external <- -1e6
    rates <- getRates(params)
    n_other <- params@initial_n_other
    dt <- 0.1

    production <- sum(getDetritusProduction(params))
    expect_lt(production, 0) # confirm the scenario is genuinely negative

    consumption <- detritus_consumption(params, n = params@initial_n, rates = rates)
    result <- detritus_dynamics(params,
        n = params@initial_n, n_other = n_other, rates = rates, dt = dt
    )

    expected_floored <- n_other$detritus * exp(-consumption * dt)
    expect_equal(result, expected_floored)
    expect_gte(result, 0)
})
