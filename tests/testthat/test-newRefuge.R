test_that("newRefuge runs without error", {
    params <- newReefParams()
    expect_error(newRefuge(params), NA)
})
