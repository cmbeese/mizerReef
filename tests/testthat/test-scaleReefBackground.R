test_that("scaleReefBackground runs without error", {
    params <- newReefParams()
    expect_error(scaleReefBackground(params, 1), NA)
})
