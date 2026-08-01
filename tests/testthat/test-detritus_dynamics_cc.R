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

test_that("detritus_dynamics_cc reduces to the no-consumption exponential decay branch when consumption is zero", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params <- setURcapacity(params)
    params@other_params$detritus$rho[] <- 0
    n <- params@initial_n
    n_other <- params@initial_n_other
    rates <- getRates(params)

    P <- sum(getDetritusProduction(params, n, rates))
    K <- params@other_params$detritus$capacity
    B0 <- detritus_biomass(params)
    expected <- B0 * exp(-1 / K * P)

    expect_equal(detritus_dynamics_cc(params, n, n_other, rates, dt = 1), expected)
})
