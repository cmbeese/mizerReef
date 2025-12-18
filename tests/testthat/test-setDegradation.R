test_that("setDegradation runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    deg_scale <- matrix(1, nrow = 2, ncol = 2)
    expect_error(setDegradation(params, deg_scale = deg_scale), NA)
})
