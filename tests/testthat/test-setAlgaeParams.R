test_that("setAlgaeParams runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(setAlgaeParams(params), NA)
})
