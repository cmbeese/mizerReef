test_that("matchReefGrowth runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(matchReefGrowth(params), NA)
})
