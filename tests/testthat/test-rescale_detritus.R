test_that("rescale_detritus runs without error", {
    params <- newReefParams()
    expect_error(rescale_detritus(params, 1), NA)
})
