test_that("getProductivity runs without error", {
    params <- newReefParams()
    expect_error(getProductivity(params), NA)
})
