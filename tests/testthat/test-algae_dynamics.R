test_that("algae_dynamics runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(algae_dynamics(params,
        n = params@initial_n,
        n_other = params@initial_n_other,
        rates = getRates(params),
        dt = 1), NA)
})

test_that("algae_dynamics floors negative production at zero", {
    # Same fix and rationale as detritus_dynamics()'s equivalent test (see
    # Finding 1 in that file): algae_growth is not currently expected to go
    # negative in practice, but algae_dynamics() shares the identical
    # convex-combination fragility, so the floor is applied defensively.
    # This directly exercises that floor via a synthetic negative growth
    # rate.
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$algae$growth <- -1e6
    rates <- getRates(params)
    n_other <- params@initial_n_other
    dt <- 0.1

    production <- sum(getAlgaeProduction(params))
    expect_lt(production, 0) # confirm the scenario is genuinely negative

    consumption <- algae_consumption(params, n = params@initial_n, rates = rates)
    result <- algae_dynamics(params,
        n = params@initial_n, n_other = n_other, rates = rates, dt = dt
    )

    expected_floored <- n_other$algae * exp(-consumption * dt)
    expect_equal(result, expected_floored)
    expect_gte(result, 0)
})
