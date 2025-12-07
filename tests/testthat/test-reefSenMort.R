test_that("reefSenMort runs without error", {
    params <- newReefParams()
    expect_error(reefSenMort(params), NA)
})
