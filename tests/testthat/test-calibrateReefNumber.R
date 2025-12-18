test_that("calibrateReefNumber runs without error", {
    data(caribbean_3_model)
    expect_error(calibrateReefNumber(caribbean_3_model), NA)
})
