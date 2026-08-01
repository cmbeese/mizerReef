test_that("algae_dynamics_cc runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(algae_dynamics_cc(params,
        n = params@initial_n,
        n_other = params@initial_n_other,
        rates = getRates(params),
        dt = 1), NA)
})

test_that("algae_dynamics_cc relaxes toward carrying capacity, not zero, when consumption is zero but production isn't", {
    # Regression test for a real bug (same pattern as detritus_dynamics_cc()):
    # the consumption == 0 branch used to return only B0 * exp(-dt/K * P)
    # (decay toward zero), dropping the `+ K * (1 - et)` relaxation-to-
    # capacity term that the general formula's own c -> 0 limit requires
    # (frac -> K, not 0). Verified independently via tiny-step Euler
    # integration of the documented ODE dB/dt = P*(1 - B/K) - c*B with c
    # fixed at 0, rather than re-deriving the closed form a second time.
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$algae$rho[] <- 0
    n <- params@initial_n
    n_other <- params@initial_n_other
    rates <- getRates(params)

    K <- 50
    params@other_params$algae$capacity <- K
    P <- sum(getAlgaeProduction(params, 0))
    B0 <- algae_biomass(params)

    dt_total <- 1
    n_steps <- 2e5
    h <- dt_total / n_steps
    B <- B0
    for (i in seq_len(n_steps)) {
        B <- B + h * (P * (1 - B / K))
    }

    expect_equal(
        algae_dynamics_cc(params, n, n_other, rates, dt = dt_total, t = 0),
        B,
        tolerance = 1e-4
    )
})

test_that("algae_dynamics_cc uses a boosted carrying capacity at a boosted timestep", {
    data(caribbean_3_model)
    data(rubble_scale)
    params <- caribbean_3_model
    params@other_dynamics$algae <- "algae_dynamics_cc"
    params <- setDegradation(params,
        trajectory = "rubble", deg_scale = rubble_scale,
        bleach_time = 2, degrade = TRUE, algae_boost = TRUE,
        algae_growth_boost = c(1), algae_capacity_boost = c(10)
    )
    # Use a small, test-fixture algae_growth rather than the bundled
    # model's real (literature-informed, much larger) production rate: with
    # capacity = 1 and dt = 1, et = exp(-dt/ka * production) underflows to
    # exactly 0 for both the boosted and unboosted ka at large production,
    # making the two next-step biomasses indistinguishable for reasons
    # unrelated to what this test checks.
    params@other_params$algae$growth <- 1
    # Zero out species abundance so algae_consumption() is 0 and
    # algae_dynamics_cc() takes its no-consumption branch
    # (et = exp(-dt/ka * production)), where ka's effect is directly
    # visible. With the model's real (very high mass-specific) consumption
    # rate, consumption dominates the dynamics regardless of ka, making the
    # capacity boost's effect unobservable in next-step biomass -- this
    # isolates the capacity term to confirm the boost is actually applied.
    n <- params@initial_n * 0
    n_other <- params@initial_n_other
    n_other$algae <- 0.5
    rates <- list()

    next_biomass_no_boost <- algae_dynamics_cc(params, n, n_other, rates,
        dt = 1, t = 0
    )
    next_biomass_boosted <- algae_dynamics_cc(params, n, n_other, rates,
        dt = 1, t = 2
    )
    expect_false(isTRUE(all.equal(next_biomass_no_boost, next_biomass_boosted)))
})
