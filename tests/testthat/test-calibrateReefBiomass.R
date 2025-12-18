test_that("calibrateReefBiomass runs without error", {
    data(caribbean_3_model)
    expect_error(calibrateReefBiomass(caribbean_3_model), NA)
})
