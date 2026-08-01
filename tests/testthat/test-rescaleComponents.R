test_that("rescaleComponents matches applying rescale_detritus then rescale_algae directly", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expected <- rescale_algae(rescale_detritus(params, 3), 2)

    result <- rescaleComponents(params, algae_factor = 2, detritus_factor = 3)
    expect_equal(result, expected)
})

test_that("rescaleComponents defaults both factors to 1 (no-op)", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    expect_equal(rescaleComponents(params), params)
})

test_that("rescaleComponents rescales algae and detritus independently", {
    data(caribbean_3_model)
    params <- caribbean_3_model
    old_algae_biomass <- algae_biomass(params)
    old_detritus_biomass <- detritus_biomass(params)

    result <- rescaleComponents(params, algae_factor = 2, detritus_factor = 1)
    expect_equal(algae_biomass(result), old_algae_biomass * 2)
    expect_equal(detritus_biomass(result), old_detritus_biomass)
})
