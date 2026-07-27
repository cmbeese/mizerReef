test_that("tuneUR_cc sets algae_growth to aout / (1 - ba/ka)", {
    # tuneUR_cc()'s documented formula: the steady-state condition for
    # algae_dynamics_cc()'s ODE dB/dt = P*(1-B/K) - c*B rearranges to
    # P = c*B / (1-B/K); since getAlgaeConsumption() already returns the
    # total (biomass-multiplied) rate c*B, this is aout / (1 - ba/ka),
    # computed here independently of tuneUR_cc()'s own call to the formula.
    data(caribbean_3_model)
    cc_params <- setURcapacity(caribbean_3_model, cap = 1.5)

    ba <- algae_biomass(cc_params)
    ka <- cc_params@other_params$algae_params$algae_capacity
    aout <- sum(getAlgaeConsumption(cc_params))
    expected_growth <- aout / (1 - ba / ka)

    result <- suppressWarnings(tuneUR_cc(cc_params))
    expect_equal(result@other_params$algae_params$algae_growth, expected_growth)
})

test_that("tuneUR_cc sets detritus external flux to dout / (1 - bd/kd) - din", {
    data(caribbean_3_model)
    cc_params <- setURcapacity(caribbean_3_model, cap = 1.5)

    bd <- detritus_biomass(cc_params)
    kd <- cc_params@other_params$detritus_params$detritus_capacity
    zeroed <- cc_params
    zeroed@other_params$detritus_params$external <- 0
    din <- sum(getDetritusProduction(zeroed))
    dout <- sum(getDetritusConsumption(cc_params))
    expected_external <- (dout / (1 - bd / kd)) - din

    result <- suppressWarnings(tuneUR_cc(cc_params))
    expect_equal(result@other_params$detritus_params$external, expected_external)
})

test_that("tuneUR_cc reaches a genuine steady state: dB/dt = 0 for both resources", {
    # Plug the tuned result back into algae_dynamics_cc()/detritus_dynamics_cc()'s
    # own documented ODE, dB/dt = P*(1-B/K) - c*B, and confirm it is exactly zero.
    # This is the check that originally caught a real bug: the old code computed
    # algae_growth <- (aout * ba) / (1 - ba/ka) and
    # external <- ((dout * bd) / (1 - bd/kd)) - din, i.e. an extra spurious
    # factor of ba/bd, since aout/dout are already biomass-multiplied totals.
    # That gave dB/dt far from zero; the fix (dividing, not multiplying, by
    # aout/dout directly) reaches exact steady state.
    data(caribbean_3_model)
    cc_params <- setURcapacity(caribbean_3_model, cap = 1.5)
    result <- suppressWarnings(tuneUR_cc(cc_params))

    ba <- algae_biomass(result)
    ka <- result@other_params$algae_params$algae_capacity
    c_A <- algae_consumption(result, n = result@initial_n, rates = getRates(result))
    P_A <- sum(getAlgaeProduction(result))
    expect_equal(P_A * (1 - ba / ka) - c_A * ba, 0)

    bd <- detritus_biomass(result)
    kd <- result@other_params$detritus_params$detritus_capacity
    c_D <- detritus_consumption(result, n = result@initial_n, rates = getRates(result))
    P_D <- sum(getDetritusProduction(result))
    expect_equal(P_D * (1 - bd / kd) - c_D * bd, 0)
})

test_that("tuneUR_cc warns when the resulting external detritus flux is negative", {
    data(caribbean_3_model)
    cc_params <- setURcapacity(caribbean_3_model, cap = 1.5)

    zeroed <- cc_params
    zeroed@other_params$detritus_params$external <- 0
    din <- sum(getDetritusProduction(zeroed))
    dout <- sum(getDetritusConsumption(cc_params))
    expect_gt(din, dout)

    expect_warning(tuneUR_cc(cc_params), "flux of external detritus is negative")
})
