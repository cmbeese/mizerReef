test_that("getSenMort runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(getSenMort(params), NA)
})
