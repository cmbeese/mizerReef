test_that("setURcapacity sets capacities to cap times the current steady-state biomass", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    ba <- algae_biomass(params)
    bd <- detritus_biomass(params)

    result <- setURcapacity(params, cap = 2)

    expect_equal(result@other_params$algae$capacity, 2 * ba)
    expect_equal(result@other_params$detritus$capacity, 2 * bd)
})

test_that("setURcapacity defaults to cap = 1.5", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    ba <- algae_biomass(params)

    result <- setURcapacity(params)
    expect_equal(result@other_params$algae$capacity, 1.5 * ba)
})

test_that("setURcapacity switches on carrying-capacity dynamics and the use_UR_cc flag", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_false(isTRUE(params@other_params$use_UR_cc))

    result <- setURcapacity(params)
    expect_true(result@other_params$use_UR_cc)
    expect_equal(result@other_dynamics$algae, "algae_dynamics_cc")
    expect_equal(result@other_dynamics$detritus, "detritus_dynamics_cc")
})
