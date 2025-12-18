test_that("setDetritusParams runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(setDetritusParams(params), NA)
})
