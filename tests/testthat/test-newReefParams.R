test_that("newReefParams returns correct class", {
    result <- newReefParams()
    expect_s4_class(result, "MizerReefParams")
})
