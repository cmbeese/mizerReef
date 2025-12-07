test_that("setDetritusParams runs without error", {
    params <- newReefParams()
    expect_error(setDetritusParams(params), NA)
})
