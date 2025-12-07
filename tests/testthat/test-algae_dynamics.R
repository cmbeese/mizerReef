test_that("algae_dynamics runs without error", {
    params <- newReefParams()
    expect_error(algae_dynamics(params, n = NULL, n_other = NULL, rates = NULL, dt = 1), NA)
})
