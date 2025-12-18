test_that("getRefuge runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(getRefuge(params), NA)
})
