test_that("plotAlgaeConsumption runs without error", {
    params <- newReefParams()
    expect_error(plotAlgaeConsumption(params), NA)
})
