test_that("getDetritusProduction runs without error", {
    params <- newReefParams()
    expect_error(getDetritusProduction(params), NA)
})
