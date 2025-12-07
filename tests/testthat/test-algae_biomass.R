test_that("algae_biomass runs without error", {
    params <- newReefParams()
    expect_error(algae_biomass(params), NA)
})
