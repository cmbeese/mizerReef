test_that("scaleReefModel runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(scaleReefModel(params, 1), NA)
})
