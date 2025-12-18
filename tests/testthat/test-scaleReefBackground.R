test_that("scaleReefBackground runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(scaleReefBackground(params, 1), NA)
})
