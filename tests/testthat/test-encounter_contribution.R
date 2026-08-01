test_that("encounter_contribution runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(encounter_contribution(params, n_other = params@initial_n_other, component = "algae"), NA)
})

test_that("encounter_contribution runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(encounter_contribution(params, n_other = params@initial_n_other, component = "detritus"), NA)
})