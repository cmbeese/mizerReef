test_that("reefSteady runs without error", {
    params <- newReefParams()
    expect_error(reefSteady(params), NA)
})
