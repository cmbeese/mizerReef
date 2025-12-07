test_that("getRefuge runs without error", {
    params <- newReefParams()
    expect_error(getRefuge(params), NA)
})
