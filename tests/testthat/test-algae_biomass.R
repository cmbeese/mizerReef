test_that("algae_biomass runs without error", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_error(algae_biomass(params), NA)
})
