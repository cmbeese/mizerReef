test_that("detritus_dynamics runs without error", {
    params <- newReefParams()
    expect_error(detritus_dynamics(params, n = NULL, n_other = NULL, rates = NULL, dt = 1), NA)
})
