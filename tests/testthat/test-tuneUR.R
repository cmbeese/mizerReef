test_that("tuneUR leaves algae_growth untouched and sets algae biomass to P_A / c_A", {
    # tuneUR() no longer retunes algae_growth to match consumption (that
    # would make algal production track grazer demand, which is not how
    # real algal primary production works). Instead it solves the fixed
    # dB_A/dt = P_A - c_A*B_A = 0 for the biomass, using algae_consumption()
    # (mass-specific, no feeding-level factor) to match exactly what
    # algae_dynamics() uses.
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_growth <- params@other_params$algae$growth
    expected_biomass <- sum(getAlgaeProduction(params)) /
        algae_consumption(params, n = params@initial_n, rates = getRates(params))

    result <- suppressWarnings(tuneUR(params))
    expect_equal(result@other_params$algae$growth, old_growth)
    expect_equal(algae_biomass(result), expected_biomass)
})

test_that("tuneUR warns and leaves algae biomass unchanged when algae consumption is zero", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@other_params$algae$rho[] <- 0
    old_biomass <- algae_biomass(params)

    expect_warning(tuneUR(params), "no finite")
    result <- suppressWarnings(tuneUR(params))
    expect_equal(algae_biomass(result), old_biomass)
})

test_that("tuneUR sets detritus external flux to consumption minus production-with-external-zeroed", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    zeroed <- params
    zeroed@other_params$detritus$external <- 0
    din <- sum(getDetritusProduction(zeroed))
    dout <- sum(getDetritusConsumption(params))
    expected_external <- dout - din

    result <- suppressWarnings(tuneUR(params))
    expect_equal(result@other_params$detritus$external, expected_external)
})

test_that("tuneUR reaches a genuine steady state: production equals consumption for both resources", {
    # A stronger check than re-deriving the formula: after tuning, plug the
    # result back into algae_dynamics()/detritus_dynamics()'s own documented
    # ODE (dB/dt = P - c*B) and confirm it lands on exactly zero.
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- suppressWarnings(tuneUR(params))

    P_A <- sum(getAlgaeProduction(result))
    c_A <- algae_consumption(result, n = result@initial_n, rates = getRates(result))
    expect_equal(P_A - c_A * algae_biomass(result), 0)

    P_D <- sum(getDetritusProduction(result))
    c_D <- detritus_consumption(result, n = result@initial_n, rates = getRates(result))
    expect_equal(P_D - c_D * detritus_biomass(result), 0)
})

test_that("tuneUR warns when the resulting external detritus flux is negative", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    zeroed <- params
    zeroed@other_params$detritus$external <- 0
    din <- sum(getDetritusProduction(zeroed))
    dout <- sum(getDetritusConsumption(params))
    expect_gt(din, dout) # caribbean_3_model triggers the negative-flux case

    expect_warning(tuneUR(params), "flux of external detritus is negative")
})
