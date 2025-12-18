test_that("detritus_biomass runs without error", {
    data(caribbean_3_model)
    expect_error(detritus_biomass(caribbean_3_model), NA)
})
