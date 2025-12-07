test_that("constant_dynamics runs without error", {
    params <- newReefParams()
    expect_error(constant_dynamics(params, n_other = NULL, component = NULL), NA)
})
