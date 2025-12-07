test_that("encounter_contribution runs without error", {
    params <- newReefParams()
    expect_error(encounter_contribution(params, n_other = NULL, component = NULL), NA)
})
