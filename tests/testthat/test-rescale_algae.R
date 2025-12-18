test_that("rescale_algae runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(rescale_algae(params, 1), NA)
})
