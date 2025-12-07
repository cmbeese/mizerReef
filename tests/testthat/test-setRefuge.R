test_that("setRefuge runs without error", {
    params <- newReefParams()
    expect_error(setRefuge(params, method = "default"), NA)
})
