test_that("rescale_algae runs without error", {
    params <- newReefParams()
    expect_error(rescale_algae(params, 1), NA)
})
