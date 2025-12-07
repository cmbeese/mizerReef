test_that("plotDetritusProduction runs without error", {
    params <- newReefParams()
    expect_error(plotDetritusProduction(params), NA)
})
