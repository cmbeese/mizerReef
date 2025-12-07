test_that("scaleReefAbundance runs without error", {
    params <- newReefParams()
    expect_error(scaleReefAbundance(params, 1), NA)
})
