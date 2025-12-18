test_that("getProductivity runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(getProductivity(params), NA)
})
