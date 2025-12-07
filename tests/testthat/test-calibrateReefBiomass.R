test_that("calibrateReefBiomass runs without error", {
    params <- newReefParams()
    expect_error(calibrateReefBiomass(params), NA)
})
