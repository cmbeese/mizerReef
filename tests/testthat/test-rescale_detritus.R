test_that("rescale_detritus runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(rescale_detritus(params, 1), NA)
})
