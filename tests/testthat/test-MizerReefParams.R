test_that("MizerReefParams returns correct class", {
    result <- MizerReefParams()
    expect_s4_class(result, "MizerReefParams")
})
