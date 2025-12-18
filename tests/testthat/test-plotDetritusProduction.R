test_that("plotDetritusProduction runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(plotDetritusProduction(params), NA)
})
