test_that("getSenMort runs without error", {
    params <- newReefParams()
    expect_error(getSenMort(params), NA)
})
