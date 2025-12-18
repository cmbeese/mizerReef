test_that("setRefuge runs without error with binned method and data frame params", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(
        start_L = 1,
        end_L = 10,
        prop_protect = 0.5
    )
    expect_error(setRefuge(params, method = "binned", method_params = method_params), NA)
})

test_that("setRefuge runs without error with binned method and list params", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- list(
        start_L = 1,
        end_L = 10,
        prop_protect = 0.5
    )
    expect_error(setRefuge(params, method = "binned", method_params = method_params), NA)
})