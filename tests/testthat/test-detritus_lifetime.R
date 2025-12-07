test_that("detritus_lifetime runs without error", {
    params <- newReefParams()
    expect_error(detritus_lifetime(params), NA)
})
