test_that("algae_dynamics_cc runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(algae_dynamics_cc(params,
        n = params@initial_n,
        n_other = params@initial_n_other,
        rates = getRates(params),
        dt = 1), NA)
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
