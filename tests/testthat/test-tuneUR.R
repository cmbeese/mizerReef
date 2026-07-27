test_that("tuneUR sets algae_growth to total algae consumption", {
    # tuneUR()'s documented formula is algae_growth <- sum(getAlgaeConsumption(params)),
    # computed here independently of tuneUR()'s own call to that formula.
    data(caribbean_3_model)
    params <- caribbean_3_model
    expected_growth <- sum(getAlgaeConsumption(params))

    result <- suppressWarnings(tuneUR(params))
    expect_equal(result@other_params$algae_params$algae_growth, expected_growth)
})

test_that("tuneUR sets detritus external flux to consumption minus production-with-external-zeroed", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    zeroed <- params
    zeroed@other_params$detritus_params$external <- 0
    din <- sum(getDetritusProduction(zeroed))
    dout <- sum(getDetritusConsumption(params))
    expected_external <- dout - din

    result <- suppressWarnings(tuneUR(params))
    expect_equal(result@other_params$detritus_params$external, expected_external)
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
    zeroed@other_params$detritus_params$external <- 0
    din <- sum(getDetritusProduction(zeroed))
    dout <- sum(getDetritusConsumption(params))
    expect_gt(din, dout) # caribbean_3_model triggers the negative-flux case

    expect_warning(tuneUR(params), "flux of external detritus is negative")
})
