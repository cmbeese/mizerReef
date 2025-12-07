test_that("setExtMortParams runs without error", {
    params <- newReefParams()
    expect_error(setExtMortParams(params), NA)
})
