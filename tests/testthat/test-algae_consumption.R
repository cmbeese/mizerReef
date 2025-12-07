test_that("algae_consumption runs without error", {
    params <- newReefParams()
    expect_error(algae_consumption(params), NA)
})
