test_that("scaleReefAbundance runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(scaleReefAbundance(params, 1), NA)
})
