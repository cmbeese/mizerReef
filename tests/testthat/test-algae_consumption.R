test_that("algae_consumption runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(algae_consumption(params), NA)
})
