test_that("setExtMortParams runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(setExtMortParams(params), NA)
})
