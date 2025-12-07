test_that("detritus_consumption runs without error", {
    params <- newReefParams()
    expect_error(detritus_consumption(params), NA)
})
