test_that("rescaleComponents runs without error", {
    params <- newReefParams()
    expect_error(rescaleComponents(params), NA)
})
