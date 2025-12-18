test_that("plotAlgaeConsumption runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(plotAlgaeConsumption(params), NA)
})
