test_that("detritus_dynamics_cc's analytic solution matches numerical integration of its own documented ODE", {
    # A real bug in the closed-form solution (wrong sign, wrong exponent,
    # etc.) would make it stop solving dB/dt = P*(1-B/K) - c*B, which this
    # would catch: independently integrate the SAME ODE with many small
    # Euler steps and confirm the two agree, rather than re-deriving the
    # same closed form a second time.
    data(caribbean_3_model)
    params <- caribbean_3_model
    n <- params@initial_n
    rates <- getRates(params)

    # The bundled model's real detritus consumption rate is astronomically
    # large (a symptom of the known, separately-tracked absolute-biomass-
    # scale calibration gap -- see algae_biomass()'s docs), which makes the
    # ODE far too numerically stiff for explicit-Euler cross-checking at
    # any feasible step count. Rescale rho down to a numerically tractable
    # magnitude first -- this only changes the size of one input, not the
    # function under test's logic, so it's still exercising the real
    # formula on real (rescaled) data.
    c_D_raw <- detritus_consumption(params, n, rates)
    params@other_params$detritus$rho <-
        params@other_params$detritus$rho * (2 / c_D_raw)

    # detritus_dynamics_cc() computes production via getDetritusProduction(params)
    # with no n/rates arguments (i.e. it lets that function recompute its own
    # fresh getRates() internally, rather than reusing the rates passed in
    # here) -- matched exactly so this test's "expected" P is computed the
    # same way the function under test actually computes it.
    P <- sum(getDetritusProduction(params))
    c_D <- detritus_consumption(params, n, rates)
    B0 <- detritus_biomass(params)
    # Likewise use a sane, well-conditioned capacity (generously larger
    # than P and c_D) rather than the real, numerically tiny,
    # biomass-derived setURcapacity() value, to keep the ODE's rate
    # constant (P/K + c_D) small enough for explicit Euler to be accurate.
    K <- 50
    params@other_params$detritus$capacity <- K
    n_other <- params@initial_n_other
    n_other$detritus <- B0

    dt_total <- 1
    rate_const <- P / K + c_D
    n_steps <- min(max(2e5, ceiling(1000 * rate_const * dt_total)), 5e6)
    h <- dt_total / n_steps
    B <- B0
    for (i in seq_len(n_steps)) {
        B <- B + h * (P * (1 - B / K) - c_D * B)
    }

    analytic <- detritus_dynamics_cc(params, n, n_other, rates, dt = dt_total)
    expect_equal(analytic, B, tolerance = 1e-3)
})

test_that("detritus_dynamics_cc relaxes toward carrying capacity, not zero, when consumption is zero but production isn't", {
    # Regression test for a real bug: the consumption == 0 branch used to
    # return only B0 * exp(-dt/K * P) (decay toward zero), dropping the
    # `+ K * (1 - et)` relaxation-to-capacity term that the general formula's
    # own c -> 0 limit requires (frac -> K, not 0). Verified independently
    # here via tiny-step Euler integration of the documented ODE
    # dB/dt = P*(1 - B/K) - c*B with c fixed at 0, rather than re-deriving
    # the closed form a second time.
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$detritus$rho[] <- 0
    n <- params@initial_n
    n_other <- params@initial_n_other
    rates <- getRates(params)

    # As in the test above, use a sane synthetic capacity rather than the
    # real, numerically tiny, biomass-derived setURcapacity() value -- the
    # real value makes dt/K * P astronomically large, which both makes
    # explicit Euler diverge AND (irrelevantly to what this test checks)
    # underflows et to exactly 0 either way.
    K <- 50
    params@other_params$detritus$capacity <- K
    P <- sum(getDetritusProduction(params))
    B0 <- detritus_biomass(params)

    dt_total <- 1
    n_steps <- 2e5
    h <- dt_total / n_steps
    B <- B0
    for (i in seq_len(n_steps)) {
        B <- B + h * (P * (1 - B / K))
    }

    expect_equal(detritus_dynamics_cc(params, n, n_other, rates, dt = dt_total),
        B,
        tolerance = 1e-4
    )
})

test_that("detritus_dynamics_cc's no-consumption branch leaves biomass unchanged when production is also zero", {
    # The one case the special branch genuinely needs to handle (avoiding
    # 0/0 in the general formula's `frac`): both production and consumption
    # zero means dB/dt = 0, so biomass should be exactly unchanged.
    # Note: detritus_dynamics_cc() computes `production` internally via
    # getDetritusProduction(params) with NO n/rates arguments (its own fresh
    # getRates() call on params' real, unzeroed abundances) rather than
    # reusing the n/rates passed in here -- so zeroing consumption's `n`
    # alone does not zero production's feces term. alpha = 1 (no waste)
    # zeroes feces regardless of abundance; sen_decomp/ext_decomp/external
    # zero out the rest.
    data(caribbean_3_model)
    params <- caribbean_3_model
    params <- setURcapacity(params)
    params@other_params$detritus$rho[] <- 0
    params@other_params$detritus$external <- 0
    params@other_params$detritus$sen_decomp <- 0
    params@other_params$detritus$ext_decomp <- 0
    params@species_params$alpha <- 1
    n <- params@initial_n * 0
    n_other <- params@initial_n_other
    rates <- getRates(params)

    expect_equal(sum(getDetritusProduction(params)), 0)
    expect_equal(
        detritus_dynamics_cc(params, n, n_other, rates, dt = 1),
        n_other$detritus
    )
})
