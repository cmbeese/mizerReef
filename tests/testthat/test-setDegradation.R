test_that("setDegradation runs without error", {
    params <- newReefParams()
    deg_scale <- matrix(1, nrow = 2, ncol = 2)
    expect_error(setDegradation(params, deg_scale = deg_scale), NA)
})
