test_that("setURcapacity runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(setURcapacity(params), NA)
})
