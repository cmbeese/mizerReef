test_that("tuneUR runs without error", {
    params <- newReefParams()
    expect_error(tuneUR(params), NA)
})
