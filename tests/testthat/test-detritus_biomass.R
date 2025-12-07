test_that("detritus_biomass runs without error", {
    params <- newReefParams()
    expect_error(detritus_biomass(params), NA)
})
