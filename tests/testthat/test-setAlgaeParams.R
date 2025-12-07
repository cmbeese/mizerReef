test_that("setAlgaeParams runs without error", {
    params <- newReefParams()
    expect_error(setAlgaeParams(params), NA)
})
