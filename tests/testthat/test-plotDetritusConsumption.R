test_that("plotDetritusConsumption runs without error", {
    params <- newReefParams()
    expect_error(plotDetritusConsumption(params), NA)
})
