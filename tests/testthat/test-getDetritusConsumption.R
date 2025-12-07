test_that("getDetritusConsumption runs without error", {
    params <- newReefParams()
    expect_error(getDetritusConsumption(params), NA)
})
