test_that("tuneUR runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(tuneUR(params), NA)
})
