test_that("calibrateReefNumber runs without error", {
    params <- newReefParams()
    expect_error(calibrateReefNumber(params), NA)
})
