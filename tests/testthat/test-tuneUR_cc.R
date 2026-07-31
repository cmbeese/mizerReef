test_that("tuneUR_cc leaves algae_growth untouched and sets algae biomass to (K*P) / (P + K*c)", {
    # tuneUR_cc() no longer retunes algae_growth to match consumption.
    # Instead it solves the fixed dB_A/dt = P*(1-B/K) - c*B = 0 for the
    # biomass, using algae_consumption() to match algae_dynamics_cc().
    data(caribbean_3_model)
    cc_params <- setURcapacity(caribbean_3_model, cap = 1.5)

    old_growth <- cc_params@other_params$algae_params$algae_growth
    ka <- cc_params@other_params$algae_params$algae_capacity
    pa <- sum(getAlgaeProduction(cc_params))
    ca <- algae_consumption(cc_params, n = cc_params@initial_n, rates = getRates(cc_params))
    expected_biomass <- (ka * pa) / (pa + ka * ca)

    result <- suppressWarnings(tuneUR_cc(cc_params))
    expect_equal(result@other_params$algae_params$algae_growth, old_growth)
    expect_equal(algae_biomass(result), expected_biomass)
})

test_that("tuneUR_cc sets detritus external flux to dout / (1 - bd/kd) - din, using the post-algae-tuning state", {
    # tuneUR_cc() tunes algae biomass first, then computes detritus's
    # production/consumption from that already-updated state (detritus
    # genuinely depends on the current system state, including algae, via
    # fish feeding level) -- so the expected values here must be computed
    # from the same post-algae-tuning intermediate state that tuneUR_cc()
    # itself uses, not from the original input params.
    data(caribbean_3_model)
    cc_params <- setURcapacity(caribbean_3_model, cap = 1.5)

    # Replicate tuneUR_cc()'s algae step to get the same intermediate state.
    ka <- cc_params@other_params$algae_params$algae_capacity
    pa <- sum(getAlgaeProduction(cc_params))
    ca <- algae_consumption(cc_params, n = cc_params@initial_n, rates = getRates(cc_params))
    algae_tuned <- cc_params
    algae_tuned@initial_n_other$algae <- (ka * pa) / (pa + ka * ca)

    bd <- detritus_biomass(algae_tuned)
    kd <- algae_tuned@other_params$detritus_params$detritus_capacity
    zeroed <- algae_tuned
    zeroed@other_params$detritus_params$external <- 0
    din <- sum(getDetritusProduction(zeroed))
    dout <- sum(getDetritusConsumption(algae_tuned))
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
