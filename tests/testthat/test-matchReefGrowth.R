test_that("matchReefGrowth runs without error", {
    params <- newReefParams()
    expect_error(matchReefGrowth(params), NA)
})
