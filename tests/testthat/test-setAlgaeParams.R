test_that("setAlgaeParams defaults match its documented values", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- setAlgaeParams(params)

    expect_equal(result@other_params$algae$growth, 2e3)
    expect_equal(result@other_params$algae$capacity, 1)
    expect_false(result@other_params$use_UR_cc)
})

test_that("setAlgaeParams stores the supplied growth and capacity values", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    result <- setAlgaeParams(params, algae_growth_initial = 42, algae_capacity = 7)

    expect_equal(result@other_params$algae$growth, 42)
    expect_equal(result@other_params$algae$capacity, 7)
})

test_that("setAlgaeParams rejects a negative growth rate or capacity", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(setAlgaeParams(params, algae_growth_initial = -1), "non-negative")
    expect_error(setAlgaeParams(params, algae_capacity = -1), "non-negative")
})

test_that("setAlgaeParams warns and zeroes interaction_algae when the column is missing", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    params@species_params$interaction_algae <- NULL

    expect_warning(result <- setAlgaeParams(params), "interaction_algae")
    expect_true(all(result@species_params$interaction_algae == 0))
})

test_that("setAlgaeParams's UR_interaction argument sets species_params columns directly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    no_sp <- nrow(params@species_params)
    custom <- rep(0.5, no_sp)

    result <- setAlgaeParams(params, UR_interaction = list(interaction_algae = custom))
    expect_equal(unname(result@species_params$interaction_algae), custom)
})
