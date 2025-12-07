test_that("algae_dynamics_cc runs without error", {
    params <- newReefParams()
    expect_error(algae_dynamics_cc(params, n = NULL, n_other = NULL, rates = NULL, dt = 1), NA)
})
