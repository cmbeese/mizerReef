test_that("getAlgaeConsumption runs without error", {
    params <- newReefParams()
    expect_error(getAlgaeConsumption(params), NA)
})
