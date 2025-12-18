test_that("detritus_consumption runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(detritus_consumption(params), NA)
})
