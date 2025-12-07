test_that("getAlgaeProduction runs without error", {
    params <- newReefParams()
    expect_error(getAlgaeProduction(params), NA)
})
