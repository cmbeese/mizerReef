test_that("reefSenMort runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(reefSenMort(params), NA)
})
