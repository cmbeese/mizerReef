test_that("newRefuge runs without error with sigmoid method", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- data.frame(L_refuge = 0.5, prop_protect = 0.3)
    expect_error(newRefuge(params, new_method = "sigmoidal", new_method_params = method_params), NA)
})

test_that("newRefuge runs without error with sigmoid method and list params", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    method_params <- list(L_refuge = 0.5, prop_protect = 0.3)
    expect_error(newRefuge(params, new_method = "sigmoidal", new_method_params = method_params), NA)
})