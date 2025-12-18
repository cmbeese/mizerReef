test_that("tuneUR_cc runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(tuneUR_cc(params), NA)
})
