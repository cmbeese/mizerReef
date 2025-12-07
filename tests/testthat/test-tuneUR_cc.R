test_that("tuneUR_cc runs without error", {
    params <- newReefParams()
    expect_error(tuneUR_cc(params), NA)
})
