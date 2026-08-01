test_that("constant_dynamics runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(constant_dynamics(params, n_other = params@initial_n_other, component = "algae"), NA)
})

test_that("constant_dynamics runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(constant_dynamics(params, n_other = params@initial_n_other, component = "detritus"), NA)
})