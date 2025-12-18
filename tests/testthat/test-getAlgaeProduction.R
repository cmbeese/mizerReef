test_that("getAlgaeProduction runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(getAlgaeProduction(params), NA)
})
